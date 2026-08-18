# hexcymatix
# // from @xoreaxeaxeax

# Overall flow:
#
#   1. Read binary file into a byte buffer.
#
#   2. process():
#      diagonal scan: for each shift i in [1, file_length), build blur[j] =
#      data[j] - data[j+i]; record runs of zeros longer than MIN_FRAGMENT as
#      Fragment(A=i+j, B=j).  The "diagonal" framing comes from the
#      self-similarity matrix: cell (r,c) = data[r] XOR data[c] is zero on the
#      main diagonal and on any parallel diagonal where a repeated byte-run
#      lives.  Iterating over shifts scans those off-diagonals one at a time.
#      Difference (subtraction) rather than equality is used so the blur values
#      carry sign information, leaving a path open to near-match detection in
#      the future (non-zero but small |blur|).  Only positive shifts are
#      scanned: i=0 is the main diagonal, where every byte trivially matches
#      itself, and shift -i is the mirror of shift +i, describing the same match
#      with A and B swapped.  A=i+j is the offset in the "shifted" copy of the
#      file; B=j is the offset in the original.  Both are absolute byte
#      positions, so each Fragment names the two locations that share content.
#
#   3. reduce_fragments():
#      two-pass merge:
#      pass 1: group by A-address, fuzzy-union overlapping B spans
#      pass 2: group by B-address, fuzzy-union overlapping A spans
#      dedup before and after (canonical key = (min(A,B), max(A,B)), so (A→B)
#      and (B→A) are treated as the same fragment).
#      The diagonal scan emits many redundant fragments: content repeated at
#      three or more offsets is reported once per pair, so the same region shows
#      up under several different shifts.  Reducing in two passes handles the
#      asymmetric case where several fragments share an A-address but span
#      slightly different B ranges (pass 1 merges those B spans), then a second
#      pass catches the converse.  The canonical (min,max) dedup key treats
#      A→B and B→A as one fragment, which matters after the merge passes have
#      pulled distinct fragments onto the same pair of addresses.  fuzz>0 lets
#      runs that are interrupted by a few differing bytes be merged into one
#      fragment, tolerating minor mutations or alignment noise in real files.
#
#   4. render():
#      circular layout: each byte offset maps to an angle on a circle; each
#      fragment is a chord (straight line) plus a perpendicular-bisector arc
#      connecting its two endpoints.  The circular layout places the start and
#      end of the file physically adjacent, which matters because many formats
#      (executables, archives) repeat header/footer structures — a linear layout
#      would make those appear at opposite extremes.  Each fragment gets two
#      marks: a chord (direct line between the two matching offsets) shows the
#      relationship at a glance, while the perpendicular-bisector arc provides a
#      second visual cue that fills the interior and curves differently
#      depending on the chord's length and angle.  Both are drawn at very low
#      alpha so that regions with many overlapping fragments accumulate
#      brightness, making dense self-similarity stand out as bright clusters
#      rather than uniform noise.

import argparse
import bisect
import io
import math
import os
import subprocess
import sys
import uuid
import webbrowser
from dataclasses import dataclass

import cairo
import numpy as np
try:
    from PIL import Image as _PILImage
    _PIL_AVAILABLE = True
except ImportError:
    _PIL_AVAILABLE = False

# Minimum run length in bytes for a match to be recorded as a fragment.
# Matches must be strictly greater than this value (i.e. 9+ bytes).
#
# TODO: this global threshold is a known limitation. The goal is variable-size
# fragments: rather than one fixed cutoff, each byte position would negotiate its
# own minimum based on what's present locally.  Sketch of the intended algorithm:
#   - collect all fragments starting at this byte; find the largest.
#   - if the largest exceeds some minimum, keep only fragments at least that size
#     (dropping short noisy ones), then advance by the *smallest* kept length.
#   - advancing by the smallest (not the largest) avoids skipping over structure
#     that shorter fragments were pointing at — fixing the boundary-alignment
#     problem where a fixed stride steps over the start of a repeated unit.
# The score/threshold that drives "continue until threshold" is still TBD.
MIN_FRAGMENT = 8

# Hack to fix flicker in invert animations where the top fragments briefly
# appear and disappear too quickly
# When --animate-invert is used, remove all fragments whose length falls within
# the N longest distinct fragment lengths (e.g. 2 removes everything at the
# two longest distinct lengths).
ANIMATE_INVERT_TRIM_TOP_LENGTHS = 10

_RENDER_SIZE = 4000
#_RENDER_SIZE = 10000

def _progress(desc, current, total, width=40):
    """Print a single-line overwriting progress bar to stdout."""
    frac = current / total if total else 1.0
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    pct = int(100 * frac)
    print(f"\r{desc}: [{bar}] {pct:3d}%", end="", flush=True)


@dataclass
class Fragment:
    A: int
    B: int
    LengthA: int
    LengthB: int


def process(data: bytes, progress: bool = False) -> list:
    """Scan all diagonal shifts to find self-similar byte runs in a binary file.

    For each shift i, computes diff[j] = data[j] - data[j+i]. Runs of
    consecutive zeros (exact byte matches) longer than MIN_FRAGMENT are
    recorded as Fragment objects. The inner loop is vectorized with NumPy.

    Only positive shifts are scanned. Shift 0 is the main diagonal, where every
    byte trivially matches itself. Negative shifts are redundant: for shift -i,
    diff[j] = data[j] - data[j-i] = -diff_pos[j-i], so the zero runs sit at
    identical positions and yield the same Fragment with A and B swapped, which
    the canonical (min, max) key treats as the same fragment.

    Args:
        data: Raw bytes of the file to analyze.

    Returns:
        List of Fragment objects representing all detected self-similar runs.
    """
    arr = np.frombuffer(data, dtype=np.uint8).astype(np.int16)
    n = len(arr)
    fragments = []
    total_shifts = max(n - 1, 0)
    step = max(1, total_shifts // 200)

    for idx, i in enumerate(range(1, n)):
        if progress and idx % step == 0:
            _progress("Identifying fragments", idx, total_shifts)

        # Compute diff only over the valid j range where 0 <= j+i < n, i.e.
        # j in [0, n-i), so the diff index is j itself.
        diff = arr[:n - i] - arr[i:]

        # Find runs of consecutive zeros via boolean edge detection.
        zeros = diff == 0
        padded = np.empty(len(zeros) + 2, dtype=bool)
        padded[0] = False
        padded[1:-1] = zeros
        padded[-1] = False
        starts  = np.where(~padded[:-1] &  padded[1:])[0]
        ends    = np.where( padded[:-1] & ~padded[1:])[0]
        lengths = ends - starts
        mask = lengths > MIN_FRAGMENT
        for s, l in zip(starts[mask], lengths[mask]):
            j = int(s)
            fragments.append(Fragment(A=i + j, B=j, LengthA=int(l), LengthB=int(l)))

    if progress:
        _progress("Identifying fragments", total_shifts, total_shifts)
        print()
    return fragments


def generate_fragments_by_range(fragments):
    """Index fragments by their canonical (min, max) offset pair.

    Using (min(A,B), max(A,B)) as the key means a match from A→B and its
    mirror B→A collapse to the same entry, so each unique repeated region
    appears only once.

    Args:
        fragments: List of Fragment objects to index.

    Returns:
        Dict mapping each (min_offset, max_offset) key to one Fragment.
    """
    by_range = {}
    for f in fragments:
        key = (min(f.A, f.B), max(f.A, f.B))
        if key not in by_range:
            by_range[key] = f
    return by_range


def remove_duplicates(fragments):
    """Remove duplicate and mirrored fragments from the list in-place.

    Two fragments are considered duplicates if they share the same canonical
    (min(A,B), max(A,B)) key, which also collapses A→B and B→A mirrors.

    Args:
        fragments: List of Fragment objects; modified in-place.
    """
    by_range = generate_fragments_by_range(fragments)
    fragments[:] = list(by_range.values())


def generate_fuzzy_union_lookup(ranges_source, fuzz):
    """Merge overlapping or nearby ranges into a minimal set of super-ranges.

    Two ranges merge when the later start is within fuzz bytes of the earlier
    end (a2 <= b1 + fuzz). Uses a sort + single linear pass, then binary
    search to map each original range to its super-range. O(n log n).

    Args:
        ranges_source: Iterable of (start, end) int tuples (inclusive).
        fuzz: Maximum gap in bytes between ranges that will still be merged.
            0 means only exactly overlapping or touching ranges merge.

    Returns:
        Dict mapping each original (start, end) tuple to the (start, end) of
        the merged super-range that absorbed it.
    """
    ranges = sorted(ranges_source)
    if not ranges:
        return {}

    merged = [list(ranges[0])]
    for a, b in ranges[1:]:
        if a <= merged[-1][1] + fuzz:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])

    merged_starts = [m[0] for m in merged]
    union_lookup = {}
    for r in ranges_source:
        idx = bisect.bisect_right(merged_starts, r[0]) - 1
        union_lookup[r] = (merged[idx][0], merged[idx][1])
    return union_lookup


def reduce_fragments(fragments, file_size, fuzz=0, progress: bool = False):
    """Collapse redundant fragments in-place using a two-pass fuzzy merge.

    The diagonal scan emits ~N fragments per real match (one per starting
    offset j), so reduction is essential, not cosmetic.

    Two passes are required because content repeated at three or more offsets
    produces fragments that overlap along either address independently.
    Pass 1: group by A-address; fuzzy-union all B spans within each group.
    Pass 2: group by B-address; fuzzy-union all A spans within each group.
    Deduplication (canonical min/max key) is applied before and after; the
    trailing pass matters most, since merging can pull distinct fragments onto
    the same pair of addresses.

    Args:
        fragments: List of Fragment objects; modified in-place.
        file_size: Total size of the source file in bytes. Currently unused
            but retained for API consistency.
        fuzz: Maximum byte gap between spans that will still be merged.
            Use fuzz > 0 to bridge matches interrupted by a few differing
            bytes, collapsing them into one fragment instead of two adjacent
            ones. Defaults to 0 (only exactly overlapping/touching spans merge).
    """
    remove_duplicates(fragments)

    by_a = {}
    for f in fragments:
        by_a.setdefault(f.A, []).append(f)
    total_a = len(by_a)
    for idx, fs in enumerate(by_a.values()):
        if progress:
            _progress("Reducing fragments", idx, total_a * 2)
        ranges_b = [(f.B, f.B + f.LengthB - 1) for f in fs]
        lookup = generate_fuzzy_union_lookup(ranges_b, fuzz)
        for f in fs:
            new_b = lookup[(f.B, f.B + f.LengthB - 1)]
            f.B = new_b[0]
            f.LengthB = new_b[1] - new_b[0] + 1

    by_b = {}
    for f in fragments:
        by_b.setdefault(f.B, []).append(f)
    for idx, fs in enumerate(by_b.values()):
        if progress:
            _progress("Reducing fragments", total_a + idx, total_a * 2)
        ranges_a = [(f.A, f.A + f.LengthA - 1) for f in fs]
        lookup = generate_fuzzy_union_lookup(ranges_a, fuzz)
        for f in fs:
            new_a = lookup[(f.A, f.A + f.LengthA - 1)]
            f.A = new_a[0]
            f.LengthA = new_a[1] - new_a[0] + 1

    if progress:
        _progress("Reducing fragments", total_a * 2, total_a * 2)
        print()
    remove_duplicates(fragments)


def fragment_brightness(fragments, file_offset, file_size, taper_window,
                        sharpness=4.0, min_brightness=0.05):
    """Compute a brightness multiplier for each fragment based on proximity to file_offset.

    Uses a Gaussian falloff exp(-sharpness * (dist/taper_window)²) floored at
    min_brightness.  Distance wraps around at the file boundaries so fragments
    near the start/end of the file taper uniformly in the circular view.

    Args:
        fragments: List of Fragment objects.
        file_offset: Current file position in bytes (the animation cursor).
        file_size: Total file size in bytes; used for circular distance.
        taper_window: Characteristic falloff width in bytes (Gaussian sigma-like).
        sharpness: Exponent scaling in the Gaussian; higher = steeper falloff.
        min_brightness: Brightness floor for distant fragments (default 0.05).

    Returns:
        numpy array of float32 with one brightness value per fragment, in [min_brightness, 1.0].
    """
    mids_a = np.array([f.A + f.LengthA / 2 for f in fragments], dtype=np.float64)
    mids_b = np.array([f.B + f.LengthB / 2 for f in fragments], dtype=np.float64)

    def _circ_dist(positions):
        d = np.abs(positions - file_offset)
        return np.minimum(d, file_size - d)

    t = np.minimum(_circ_dist(mids_a), _circ_dist(mids_b)) / taper_window
    return np.maximum(np.exp(-sharpness * t * t), min_brightness).astype(np.float32)


def _fragment_endpoints(f, file_size, cx, cy, R):
    """Compute the canvas (x, y) coordinates for both ends of a fragment.

    Each byte offset is mapped to an angle via offset/file_size * 2π, then
    projected onto the circle of radius R centered at (cx, cy).

    Args:
        f: Fragment with A and B byte offsets.
        file_size: Total file size in bytes; used to normalize offsets to [0, 1].
        cx: X coordinate of the circle center in pixels.
        cy: Y coordinate of the circle center in pixels.
        R: Radius of the circle in pixels.

    Returns:
        Tuple (ax, ay, bx, by) — canvas coordinates for the A and B endpoints.
    """
    TWO_PI = 2 * math.pi
    ax = cx + R * math.cos(TWO_PI * f.A / file_size)
    ay = cy + R * math.sin(TWO_PI * f.A / file_size)
    bx = cx + R * math.cos(TWO_PI * f.B / file_size)
    by = cy + R * math.sin(TWO_PI * f.B / file_size)
    return ax, ay, bx, by


def _draw_chord(ctx, f, file_size, cx, cy, R):
    """Stroke a straight chord connecting the two endpoints of a fragment.

    Args:
        ctx: Active cairo.Context; caller sets line width and source color.
        f: Fragment with A and B byte offsets.
        file_size: Total file size in bytes.
        cx: X coordinate of the circle center in pixels.
        cy: Y coordinate of the circle center in pixels.
        R: Radius of the circle in pixels.
    """
    ax, ay, bx, by = _fragment_endpoints(f, file_size, cx, cy, R)
    ctx.move_to(ax, ay)
    ctx.line_to(bx, by)
    ctx.stroke()


def _draw_arc(ctx, f, file_size, cx, cy, R):
    """Stroke a perpendicular-bisector arc through both endpoints of a fragment.

    Finds the arc whose circle passes through both endpoints with its center
    on the perpendicular bisector, computed by intersecting tangent lines at
    each point. Very large or near-degenerate arcs are scaled down via an
    arc_scale clamp to stay renderable. Silently skips degenerate cases.

    Args:
        ctx: Active cairo.Context; caller sets line width and source color.
        f: Fragment with A and B byte offsets.
        file_size: Total file size in bytes.
        cx: X coordinate of the circle center in pixels.
        cy: Y coordinate of the circle center in pixels.
        R: Radius of the circle in pixels.
    """
    # Arc center is on the perpendicular bisector of the chord, computed by
    # intersecting the two tangent lines at each endpoint.
    ax, ay, bx, by = _fragment_endpoints(f, file_size, cx, cy, R)
    dy_a = ay - cy
    dy_b = by - cy
    if abs(dy_a) < 1e-10 or abs(dy_b) < 1e-10:
        return
    slope_a = -(ax - cx) / dy_a
    slope_b = -(bx - cx) / dy_b
    if abs(slope_a - slope_b) < 1e-10:
        return
    intercept_a = ay - slope_a * ax
    intercept_b = by - slope_b * bx
    arc_cx = (intercept_a - intercept_b) / (slope_b - slope_a)
    arc_cy = slope_a * arc_cx + intercept_a
    arc_radius = math.sqrt((ax - arc_cx) ** 2 + (ay - arc_cy) ** 2)
    if arc_radius < 1e-10:
        return
    arc_scale = 1.0
    if arc_radius > 1_000_000:
        arc_scale = 1_000_000 / arc_radius
    arc_radius *= arc_scale
    arc_cx = cx * (1 - arc_scale) + arc_cx * arc_scale
    arc_cy = cy * (1 - arc_scale) + arc_cy * arc_scale
    ctx.arc(arc_cx, arc_cy, arc_radius, 0, 2 * math.pi)
    ctx.stroke()


def _draw_linear_edge_arc(ctx, edge_y, x1, x2, y_border, top):
    """Stroke a decorative bulge-arc at one edge of the linear map.

    Uses the same geometry as the framing arcs: the arc extends exactly
    y_border pixels past the bar edge when the span is large enough. For
    spans smaller than y_border, bulge is clamped to half_span, producing
    a semicircle (offset=0) that still extends as far as it can.

    Derivation: given half_span H and desired bulge B, the circle center
    sits at offset = (H² − B²) / (2·B) inward from the bar edge, with
    radius = sqrt(offset² + H²). Setting B=H gives offset=0 (semicircle).
    """
    dx = abs(x1 - x2)
    if dx < 1e-3:
        return
    half_span = dx / 2
    bulge = min(half_span, y_border)
    center_offset = (half_span ** 2 - bulge ** 2) / (2 * bulge)
    arc_radius = math.sqrt(center_offset ** 2 + half_span ** 2)
    if arc_radius < 1e-3:
        return
    half_angle = math.asin(min(1.0, half_span / arc_radius))
    cx = (x1 + x2) / 2
    if top:
        # Arc opens upward; center is below the bar (Y increases downward).
        ctx.arc(cx, edge_y + center_offset, arc_radius, -math.pi / 2 - half_angle, -math.pi / 2 + half_angle)
    else:
        # Arc opens downward; center is above the bar.
        ctx.arc(cx, edge_y - center_offset, arc_radius, math.pi / 2 - half_angle, math.pi / 2 + half_angle)
    ctx.stroke()


def _draw_linear_fragment(ctx, f, file_size, bar_a_y, bar_b_y, bar_left, bar_width, y_border, alpha_fill, alpha_arc):
    """Draw one fragment as two crossed quadrilaterals plus four edge arcs.

    The A-span and B-span are both mapped to X positions on the canvas.
    bar_a_y is the top bar Y, bar_b_y the bottom bar Y. The two quads ("forward"
    and "cross") connect the spans across the gap. Four decorative arcs are drawn
    at the top and bottom edges at each span boundary, extending up to
    y_border pixels past the bar (clamped to a semicircle for short spans).
    """
    a_x1 = bar_left + bar_width * f.A / file_size
    a_x2 = bar_left + bar_width * (f.A + f.LengthA) / file_size
    b_x1 = bar_left + bar_width * f.B / file_size
    b_x2 = bar_left + bar_width * (f.B + f.LengthB) / file_size

    ctx.set_source_rgba(1, 1, 1, alpha_fill)

    # Forward quad: A-start→B-start, A-end→B-end
    ctx.move_to(a_x1, bar_a_y)
    ctx.line_to(b_x1, bar_b_y)
    ctx.line_to(b_x2, bar_b_y)
    ctx.line_to(a_x2, bar_a_y)
    ctx.close_path()
    ctx.fill()

    # Cross quad: A-start→B-end, A-end→B-start
    ctx.move_to(b_x1, bar_a_y)
    ctx.line_to(a_x1, bar_b_y)
    ctx.line_to(a_x2, bar_b_y)
    ctx.line_to(b_x2, bar_a_y)
    ctx.close_path()
    ctx.fill()

    ctx.set_source_rgba(1, 1, 1, alpha_arc)
    _draw_linear_edge_arc(ctx, bar_a_y, a_x1, b_x1, y_border, top=True)
    _draw_linear_edge_arc(ctx, bar_a_y, a_x2, b_x2, y_border, top=True)
    _draw_linear_edge_arc(ctx, bar_b_y, a_x1, b_x1, y_border, top=False)
    _draw_linear_edge_arc(ctx, bar_b_y, a_x2, b_x2, y_border, top=False)


def draw_circle_fragment(ctx, f, file_size, chord_alpha, arc_alpha, cx, cy, R):
    """Draw one fragment as a chord and a perpendicular-bisector arc.

    Each byte offset maps to an angle on the circle (offset/file_size * 2π).
    The chord is drawn at chord_alpha and the arc at arc_alpha, both in white,
    so dense overlapping regions accumulate brightness.

    Args:
        ctx: Active cairo.Context.
        f: Fragment with A and B byte offsets.
        file_size: Total file size in bytes.
        chord_alpha: Alpha value (0–1) for the chord stroke.
        arc_alpha: Alpha value (0–1) for the arc stroke.
        cx: X coordinate of the circle center in pixels.
        cy: Y coordinate of the circle center in pixels.
        R: Radius of the circle in pixels.
    """
    ctx.set_source_rgba(1, 1, 1, chord_alpha)
    _draw_chord(ctx, f, file_size, cx, cy, R)
    ctx.set_source_rgba(1, 1, 1, arc_alpha)
    _draw_arc(ctx, f, file_size, cx, cy, R)


NAMED_RESOLUTIONS = {
    "480p":  (854,  480,  0.25),
    "720p":  (1280, 720,  0.25),
    "1080p": (1920, 1080, 0.25),
    "1440p": (2560, 1440, 0.25),
    "4k":    (3840, 2160, 0.25),
    "2160p": (3840, 2160, 0.25),
    "uwfhd": (2560, 1080, 0.50),
    "suwfhd": (3840, 1080, 0.50),
}


def parse_resolution(s):
    """Parse a resolution string into a (width, height, circle_fill) tuple.

    Accepts:
        - Named presets: 480p, 720p, 1080p, 1440p, 4k, 2160p, uwfhd, suwfhd
        - WxH format:    e.g. "1920x1080" or "2000x2000"

    Args:
        s: Resolution string from the command line.

    Returns:
        (width, height, circle_fill) tuple; circle_fill is the preset default
        fill fraction (0.25 for standard presets, 0.50 for ultra-wide presets,
        0.25 for raw WxH input).

    Raises:
        argparse.ArgumentTypeError on unrecognised input.
    """
    lower = s.lower()
    if lower in NAMED_RESOLUTIONS:
        return NAMED_RESOLUTIONS[lower]
    parts = lower.split("x")
    if len(parts) == 2:
        try:
            w, h = int(parts[0]), int(parts[1])
            if w > 0 and h > 0:
                return (w, h, 0.25)
        except ValueError:
            pass
    raise argparse.ArgumentTypeError(
        f"unrecognised resolution '{s}'; use WxH (e.g. 2000x2000) "
        f"or one of: {', '.join(NAMED_RESOLUTIONS)}"
    )


def positive_int(s):
    """Parse a strictly positive integer for argparse.

    Guards --top, where a value below 1 would index fragments[top_n - 1] from
    the end of the list and silently keep every fragment instead of none.

    Raises:
        argparse.ArgumentTypeError on non-integer or non-positive input.
    """
    try:
        v = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{s}' is not an integer")
    if v < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or greater, got {v}")
    return v


def render(fragments, file_size, progress: bool = False, circle_fill: float = 0.25, brightness=None, draw_arcs: bool = True):
    """Render all fragments as a circular self-similarity map.

    Always draws on a 4000×4000 black canvas. Use downscale() afterwards to
    resize to the desired output resolution. Each fragment is a chord plus a
    perpendicular-bisector arc at low alpha so overlapping fragments accumulate
    brightness. Three radial vignette passes darken the rim so the circle reads
    as a sphere and outer arcs fade to black.

    Args:
        fragments: List of Fragment objects to draw.
        file_size: Total file size in bytes; maps offsets to circle angles.
        progress: If True, print a progress bar to stdout.

    Returns:
        A cairo.ImageSurface (FORMAT_RGB24, 4000×4000) with the rendered visualization.
    """
    WIDTH, HEIGHT = _RENDER_SIZE, _RENDER_SIZE

    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    ctx.set_source_rgb(0, 0, 0)
    ctx.paint()

    cx, cy = WIDTH / 2.0, HEIGHT / 2.0
    circle_diameter = min(WIDTH, HEIGHT) * circle_fill
    R = circle_diameter / 2.0

    alpha_scale = min(1.0, math.sqrt(10000 / max(len(fragments), 1)))
    arc_alpha   = (22 / 255) * alpha_scale
    chord_alpha = (50 / 255) * alpha_scale

    ctx.set_line_width(1)
    ctx.set_source_rgba(1, 1, 1, 128 / 255)
    ctx.arc(cx, cy, R, 0, 2 * math.pi)
    ctx.stroke()

    total = len(fragments)
    step = max(1, total // 200)
    if brightness is None:
        # Two-pass rendering: set source color once per pass instead of once per
        # fragment, reducing Cairo API calls from 2*N to 2 for large fragment sets.
        ctx.set_source_rgba(1, 1, 1, chord_alpha)
        for idx, f in enumerate(fragments):
            if progress and idx % step == 0:
                _progress("Rendering", idx, total * 2)
            _draw_chord(ctx, f, file_size, cx, cy, R)

        if draw_arcs:
            ctx.set_source_rgba(1, 1, 1, arc_alpha)
            for idx, f in enumerate(fragments):
                if progress and idx % step == 0:
                    _progress("Rendering", total + idx, total * 2)
                _draw_arc(ctx, f, file_size, cx, cy, R)
    else:
        for idx, f in enumerate(fragments):
            if progress and idx % step == 0:
                _progress("Rendering", idx, total * 2)
            b = float(brightness[idx])
            ctx.set_source_rgba(1, 1, 1, chord_alpha * b)
            _draw_chord(ctx, f, file_size, cx, cy, R)

        if draw_arcs:
            for idx, f in enumerate(fragments):
                if progress and idx % step == 0:
                    _progress("Rendering", total + idx, total * 2)
                b = float(brightness[idx])
                ctx.set_source_rgba(1, 1, 1, arc_alpha * b)
                _draw_arc(ctx, f, file_size, cx, cy, R)

    if progress:
        _progress("Rendering", total * 2, total * 2)
        print()

    ctx.set_source_rgba(0, 0, 0, 250 / 255)
    ctx.set_line_width(3)
    ctx.arc(cx, cy, R, 0, 2 * math.pi)
    ctx.stroke()
    ctx.set_line_width(1)

    # Sphere vignette pass 1: transparent center → black rim, clipped to ellipse.
    # GDI+ PathGradientBrush clips to its path; we must do the same or the corners
    # of the bounding box get painted black and cover arcs outside the circle.
    ctx.save()
    ctx.arc(cx, cy, R, 0, 2 * math.pi)
    ctx.clip()
    grad = cairo.RadialGradient(cx, cy, R * 0.9, cx, cy, R)
    grad.add_color_stop_rgba(0, 0, 0, 0, 0)
    grad.add_color_stop_rgba(1, 0, 0, 0, 1)
    ctx.set_source(grad)
    ctx.paint()
    ctx.restore()

    # Sphere vignette pass 2: deeper darkening toward rim (FocusScales 0.4).
    ctx.save()
    ctx.arc(cx, cy, R, 0, 2 * math.pi)
    ctx.clip()
    grad2 = cairo.RadialGradient(cx, cy, R * 0.4, cx, cy, R)
    grad2.add_color_stop_rgba(0, 0, 0, 0, 0)
    grad2.add_color_stop_rgba(1, 0, 0, 0, 200 / 255)
    ctx.set_source(grad2)
    ctx.paint()
    ctx.restore()

    # Fade to black as approach edge of canvas
    fade_start = 0.5  # fraction of R*4 where outer fade begins (0=center, 1=edge)
    #d = R * 4
    d = min(WIDTH, HEIGHT) / 2
    grad3 = cairo.RadialGradient(cx, cy, 0, cx, cy, d)
    grad3.add_color_stop_rgba(0,                                       0, 0, 0, 0)
    grad3.add_color_stop_rgba(fade_start,                              0, 0, 0, 0)
    grad3.add_color_stop_rgba(fade_start + (1 - fade_start) * 0.4,     0, 0, 0, 0.7)
    grad3.add_color_stop_rgba(fade_start + (1 - fade_start) * 0.6,     0, 0, 0, 0.9)
    grad3.add_color_stop_rgba(fade_start + (1 - fade_start) * 0.8,     0, 0, 0, 0.95)
    grad3.add_color_stop_rgba(1.0,                                     0, 0, 0, 1.0)
    ctx.set_source(grad3)
    ctx.paint()

    return surface


def render_linear(fragments, file_size, progress: bool = False, tw: int = _RENDER_SIZE, th: int = _RENDER_SIZE, draw_arcs: bool = True):
    """Render all fragments as a linear (parallel-bars) self-similarity map.

    The file is represented as two horizontal bars (top = A-offsets, bottom =
    B-offsets). Each fragment is drawn as two crossed quadrilaterals connecting
    its spans on each bar, plus four decorative bulge-arcs at the bar edges.
    Mirrors GenerateFragmentMapLinear / DrawLinearFragment from the C# source.

    Args:
        fragments: List of Fragment objects to draw.
        file_size: Total file size in bytes; maps offsets to X positions.
        progress: If True, print a progress bar to stdout.
        tw: Target output width in pixels (used to compute the cropped layout).
        th: Target output height in pixels (used to compute the cropped layout).

    Returns:
        A cairo.ImageSurface (FORMAT_RGB24, _RENDER_SIZE×_RENDER_SIZE).
    """
    WIDTH, HEIGHT = _RENDER_SIZE, _RENDER_SIZE

    # Compute the effective visible rectangle on the square render canvas after
    # downscale()'s fill+crop to (tw, th).  The longer target dimension maps to
    # the full _RENDER_SIZE; the shorter one is centered and occupies a fraction.
    if tw >= th:
        eff_w = WIDTH
        eff_h = HEIGHT * th / tw
    else:
        eff_w = WIDTH * tw / th
        eff_h = HEIGHT

    cx = WIDTH / 2
    cy = HEIGHT / 2

    # Leave at least 5% of the visible area as a blank border on every edge.
    # usable_h / usable_w are the dimensions available after reserving borders.
    # Within usable_h the 1:1:1 split is preserved (arc : bars : arc), so
    # y_border = usable_h/3 and bars sit ±usable_h/6 from the canvas center.
    edge_border = 0.05
    usable_h = eff_h * (1 - 2 * edge_border)
    usable_w = eff_w * (1 - 2 * edge_border)
    y_border = usable_h / 3
    x_border = (eff_w - usable_w) / 2  # = eff_w * edge_border, kept symmetric

    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(0, 0, 0)
    ctx.paint()

    bar_a_y  = cy - usable_h / 6
    bar_b_y  = cy + usable_h / 6
    bar_left  = cx - eff_w / 2 + x_border
    bar_right = cx + eff_w / 2 - x_border
    bar_width = bar_right - bar_left
    bar_height = bar_b_y - bar_a_y

    alpha_scale = min(1.0, math.sqrt(10000 / max(len(fragments), 1)))
    alpha_fill  = (20 / 255) * alpha_scale
    alpha_arc   = (8 / 255) * alpha_scale if draw_arcs else 0.0

    ctx.set_line_width(1)

    total = len(fragments)
    step  = max(1, total // 200)
    for idx, f in enumerate(fragments):
        if progress and idx % step == 0:
            _progress("Rendering", idx, total)
        _draw_linear_fragment(ctx, f, file_size, bar_a_y, bar_b_y, bar_left, bar_width, y_border, alpha_fill, alpha_arc)

    if progress:
        _progress("Rendering", total, total)
        print()

    # Rectangle outline (horizontal bars + left/right connectors)
    ctx.set_source_rgba(1, 1, 1, 128 / 255)
    ctx.set_line_width(2)
    ctx.rectangle(bar_left, bar_a_y, bar_width, bar_height)
    ctx.stroke()
    ctx.set_line_width(1)

    return surface


def downscale(src, tw, th):
    """Downscale a 4000×4000 cairo surface to (tw, th) using fill + center-crop.

    Scales the source so the shorter target dimension is fully covered (fill),
    then crops the longer dimension symmetrically around the center. Uses
    Pillow's Lanczos filter for high-quality downsampling.

    Args:
        src: cairo.ImageSurface produced by render() (4000×4000).
        tw: Target width in pixels.
        th: Target height in pixels.

    Returns:
        A new cairo.ImageSurface (FORMAT_RGB24, tw×th).
    """
    if _PIL_AVAILABLE:
        return _downscale_pillow(src, tw, th)
    print("Warning: Pillow not available; falling back to cairo bilinear downscale. "
          "Install Pillow for higher quality output.", file=sys.stderr)
    return _downscale_cairo(src, tw, th)


def _downscale_pillow(src, tw, th):
    # cairo FORMAT_RGB24 stores pixels as 32-bit words (BGRX on little-endian),
    # so we reinterpret as RGBA and drop the alpha channel to get RGB.
    buf = src.get_data()
    img = _PILImage.frombuffer("RGBA", (_RENDER_SIZE, _RENDER_SIZE), bytes(buf), "raw", "BGRA", 0, 1)
    img = img.convert("RGB")

    scale = max(tw / _RENDER_SIZE, th / _RENDER_SIZE)
    scaled_w = round(_RENDER_SIZE * scale)
    scaled_h = round(_RENDER_SIZE * scale)
    img = img.resize((scaled_w, scaled_h), _PILImage.LANCZOS)

    left = (scaled_w - tw) // 2
    top  = (scaled_h - th) // 2
    img = img.crop((left, top, left + tw, top + th))

    out = cairo.ImageSurface(cairo.FORMAT_RGB24, tw, th)
    out_buf = out.get_data()
    # Convert back to BGRX for cairo.
    out_buf[:] = img.tobytes("raw", "BGRX")
    return out


def _surface_to_pil(surface):
    """Convert a cairo ImageSurface (RGB24 or ARGB32) to a PIL RGB Image."""
    w = surface.get_width()
    h = surface.get_height()
    buf = surface.get_data()
    img = _PILImage.frombuffer("RGBA", (w, h), bytes(buf), "raw", "BGRA", 0, 1)
    return img.convert("RGB")


def _save_jpg(surface, png_path):
    """Save *surface* as an 80% quality JPEG alongside *png_path*."""
    if not _PIL_AVAILABLE:
        print("Warning: Pillow not available; skipping JPEG output.", file=sys.stderr)
        return
    jpg_path = os.path.splitext(png_path)[0] + ".jpg"
    img = _surface_to_pil(surface)
    img.save(jpg_path, "JPEG", quality=80)
    print(f"Saved: {jpg_path}")


def _downscale_cairo(src, tw, th):
    scale = max(tw / _RENDER_SIZE, th / _RENDER_SIZE)
    ox = (_RENDER_SIZE * scale - tw) / 2
    oy = (_RENDER_SIZE * scale - th) / 2

    out = cairo.ImageSurface(cairo.FORMAT_RGB24, tw, th)
    ctx = cairo.Context(out)
    ctx.translate(-ox, -oy)
    ctx.scale(scale, scale)
    ctx.set_source_surface(src, 0, 0)
    ctx.get_source().set_filter(cairo.Filter.BILINEAR)
    ctx.paint()
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fragment identification and circular visualizer.")
    parser.add_argument("files", nargs="+", metavar="file", help="Binary file(s) to analyze")
    parser.add_argument(
        "--resolution",
        type=parse_resolution,
        default=(1920, 1080, 0.25),
        metavar="RES",
        help="Output resolution as WxH or a named preset (480p, 720p, 1080p, 1440p, 4k). Default: 1080p",
    )
    parser.add_argument(
        "--ultra",
        action="store_const",
        const=(4000, 4000, 0.25),
        dest="resolution",
        help="Shorthand for --resolution 4000x4000",
    )
    parser.add_argument(
        "--top",
        type=positive_int,
        nargs="+",
        default=None,
        metavar="N",
        help="Keep only the N longest fragments before rendering; repeat or space-separate for multiple renders",
    )
    parser.add_argument(
        "--animate",
        action="store_true",
        help="Render an MP4 animation (first 10%%, 20%%, ... 100%% of fragments)",
    )
    parser.add_argument(
        "--animate-2",
        action="store_true",
        help="Render an MP4 animation revealing fragments shortest-first; "
             "first frame shows just the circle, last frame shows all fragments",
    )
    parser.add_argument(
        "--animate-fade",
        action="store_true",
        help="With --animate-2: fade older fragments via Gaussian falloff so "
             "only the most recently revealed fragments are at full brightness",
    )
    parser.add_argument(
        "--animate-invert",
        action="store_true",
        help="With --animate-2: after all fragments are revealed, hold for 10%% "
             "of frames then remove fragments in reverse order back to the circle",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=10,
        metavar="N",
        help="Number of animation frames (default: 10, used with --animate)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=20,
        metavar="N",
        help="Frames per second for the output video (default: 20, used with --animate)",
    )
    parser.add_argument(
        "--taper-window",
        type=float,
        default=0.10,
        metavar="FRAC",
        help="Animation taper window as a fraction of file size (default: 0.10). "
             "Controls the characteristic Gaussian falloff width.",
    )
    parser.add_argument(
        "--taper-sharpness",
        type=float,
        default=4.0,
        metavar="K",
        help="Gaussian sharpness for animation brightness taper (default: 4.0). "
             "Higher values produce a steeper falloff.",
    )
    parser.add_argument(
        "--circular",
        action="store_true",
        help="Render as circular layout (default if neither --circular nor --linear is given)",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Render as parallel-bars layout",
    )
    parser.add_argument(
        "--circle-fill",
        type=float,
        default=None,
        metavar="FRACTION",
        help="Fraction of the canvas diameter used by the circle (e.g. 0.25). "
             "Overrides the resolution preset default.",
    )
    parser.add_argument(
        "--no-arcs",
        action="store_true",
        help="Skip drawing arcs; render chords only (faster, less visual noise)",
    )
    parser.add_argument(
        "--jpg",
        action="store_true",
        help="Also save an 80%% quality JPEG alongside each PNG output",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        metavar="DIR",
        help="Directory to write output files into (default: output)",
    )
    args = parser.parse_args()

    # Reject flag combinations that would otherwise be silently ignored.
    if args.animate and args.animate_2:
        parser.error("--animate and --animate-2 are mutually exclusive")
    if (args.animate or args.animate_2) and args.linear:
        parser.error("--animate and --animate-2 render the circular layout only; "
                     "--linear is not supported")
    if (args.animate_fade or args.animate_invert) and not args.animate_2:
        parser.error("--animate-fade and --animate-invert require --animate-2")

    tw, th, res_fill = args.resolution
    circle_fill = args.circle_fill if args.circle_fill is not None else res_fill
    # Convert circle_fill (fraction of min output dimension) to the fraction of
    # the high-res render canvas that the circle should occupy, so that after
    # downscale's fill+crop the circle lands at the requested size.
    render_fill = circle_fill * min(tw, th) / max(tw, th)

    for path in args.files:
        try:
            data = open(path, "rb").read()
        except OSError as e:
            print(f"Warning: could not open {path}: {e}")
            continue
        n = len(data)
        print(f"Processing {path} ({n} bytes)...")

        fragments = process(data, progress=True)
        print(f"Found {len(fragments)} raw fragments")

        reduce_fragments(fragments, n, fuzz=0, progress=True)
        print(f"Reduced to {len(fragments)} fragments")

        # Sort by length once so each --top value can slice without re-sorting.
        if args.top is not None:
            fragments.sort(key=lambda f: max(f.LengthA, f.LengthB), reverse=True)

        top_values = args.top if args.top is not None else [None]
        for top_n in top_values:
            if top_n is not None:
                before_top = len(fragments)
                render_frags = fragments[:before_top]  # start with all
                if before_top > top_n:
                    # Expand the cutoff to include all fragments tied at the boundary
                    # length, so --top never arbitrarily drops half a tie group.
                    cutoff_len = max(fragments[top_n - 1].LengthA, fragments[top_n - 1].LengthB)
                    render_frags = [f for f in fragments if max(f.LengthA, f.LengthB) >= cutoff_len]
                    print(f"Warning: --top truncated {before_top} fragments to {len(render_frags)} "
                          f"(cutoff length: {cutoff_len}); consider raising --top")
                print(f"Keeping top {len(render_frags)} fragments by length")
            else:
                render_frags = fragments[:]

            # Sort by A*4 + B: primarily by A, with B as a secondary key. The factor 4
            # controls the width of the B-window within which fragments at nearby A values
            # can interleave: two fragments swap order only when their B values differ by
            # more than 4 * |ΔA|. A larger factor narrows this window (A dominates more
            # strongly); a smaller factor widens it (B has more influence). At 4, only
            # fragments whose B offsets are within ~4 bytes per unit of A-distance will
            # reorder relative to each other — in practice A dominates and B breaks ties.
            render_frags.sort(key=lambda f: f.A * 4 + f.B)
            #render_frags.sort(key=lambda f: f.A)

            if args.animate:
                # Each frame sweeps a file-offset cursor across the file.  All fragments
                # are rendered every frame; brightness is 1.0 near the cursor and tapers
                # to 10% for fragments far from it.  The cursor wraps around so fragments
                # near the file boundaries get a uniform taper on both sides.
                taper_window = int(n * args.taper_window)

                os.makedirs(args.output_dir, exist_ok=True)
                base = os.path.basename(path)
                res_tag = f"_{tw}x{th}"
                top_tag = f"_top{top_n}" if top_n is not None else ""
                out_path = os.path.join(args.output_dir, f"{base}_circular{res_tag}{top_tag}_{uuid.uuid4()}.mp4")
                # yuv420p requires even dimensions
                enc_w = tw + (tw % 2)
                enc_h = th + (th % 2)
                cmd = [
                    "ffmpeg", "-y",
                    "-f", "image2pipe", "-vcodec", "png", "-r", str(args.fps), "-i", "-",
                    "-vcodec", "libx264", "-pix_fmt", "yuv420p",
                    "-vf", f"scale={enc_w}:{enc_h}",
                    out_path,
                ]
                proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)
                for i in range(args.frames):
                    file_offset = int(n * i / args.frames)
                    print(f"Rendering frame {i + 1}/{args.frames} (offset {file_offset}/{n})...")
                    bvals = fragment_brightness(render_frags, file_offset, n, taper_window,
                                               sharpness=args.taper_sharpness)
                    surface = render(render_frags, n, progress=True, circle_fill=render_fill, brightness=bvals, draw_arcs=not args.no_arcs)
                    if (tw, th) != (_RENDER_SIZE, _RENDER_SIZE):
                        surface = downscale(surface, tw, th)
                    buf = io.BytesIO()
                    surface.write_to_png(buf)
                    proc.stdin.write(buf.getvalue())
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    print(f"ffmpeg exited with code {proc.returncode}", file=sys.stderr)
                    sys.exit(1)
                print(f"Saved: {out_path}")
                webbrowser.open(f"file://{os.path.abspath(out_path)}")
            elif args.animate_2:
                # Sort fragments shortest-first; randomize within the same length
                # using numpy lexsort so ties get a stable-but-random ordering.
                lengths = np.array([max(f.LengthA, f.LengthB) for f in render_frags])
                noise   = np.random.default_rng().random(len(render_frags))
                order   = np.lexsort((noise, lengths))
                sorted_frags = [render_frags[i] for i in order]

                if args.animate_invert:
                    sorted_lengths = sorted({max(f.LengthA, f.LengthB) for f in sorted_frags}, reverse=True)
                    trim_lengths   = set(sorted_lengths[:ANIMATE_INVERT_TRIM_TOP_LENGTHS])
                    before_trim    = len(sorted_frags)
                    sorted_frags   = [f for f in sorted_frags if max(f.LengthA, f.LengthB) not in trim_lengths]
                    trimmed_count  = before_trim - len(sorted_frags)
                    if trimmed_count:
                        print(
                            f"Warning: removed {trimmed_count} fragment(s) whose length fell in the "
                            f"top {ANIMATE_INVERT_TRIM_TOP_LENGTHS} distinct lengths "
                            f"({sorted(trim_lengths, reverse=True)}) for --animate-invert",
                            file=sys.stderr,
                        )

                total_frags  = len(sorted_frags)

                taper_window_count = max(1, int(total_frags * args.taper_window))

                # Build the sequence of (k, direction) pairs where k is the
                # number of fragments to show and direction is 1=forward or
                # -1=reverse (used by --animate-fade to know which end is newest).
                #
                # --animate-invert splits args.frames across three phases:
                #   forward : args.frames frames, ramp 0→all
                #   hold    : max(1, round(args.frames * 0.2)) frames at all
                #   reverse : same count as forward, ramp all→0
                # Without --animate-invert the sequence is just the forward ramp.
                if args.animate_invert:
                    fwd_frames  = args.frames
                    #hold_frames = max(1, round(args.frames * 0.2))
                    hold_frames = 1
                    rev_frames  = args.frames
                    frame_ks = []
                    # forward ramp: 0 → total_frags
                    for i in range(fwd_frames):
                        k = int(total_frags * i / (fwd_frames - 1)) if fwd_frames > 1 else total_frags
                        frame_ks.append((k, 1))
                    # hold: oscillate back and forth at the same per-frame rate as
                    # the forward pass using a fixed number of oscillations.  The
                    # amplitude is derived so that the rate exactly matches forward:
                    #   rate = total_frags / fwd_frames  (frags per frame)
                    #   half_period = hold_frames / (2 * NUM_HOLD_OSC)
                    #   amplitude = rate * half_period
                    # This means a 10% hold with 5 oscillations traverses the last
                    # 1% of the fragment sequence per half-cycle.
                    NUM_HOLD_OSC  = 1
                    osc_amplitude = (total_frags / fwd_frames) * (hold_frames / (2 * NUM_HOLD_OSC)) if fwd_frames > 1 else 0
                    for i in range(hold_frames):
                        t     = i / hold_frames if hold_frames > 1 else 0.0
                        osc_t = (t * 2 * NUM_HOLD_OSC) % 2.0  # 0..2 per cycle
                        if osc_t < 1.0:
                            displacement = osc_amplitude * osc_t         # descending
                            direction    = -1
                        else:
                            displacement = osc_amplitude * (2.0 - osc_t) # ascending
                            direction    = 1
                        k = max(0, min(total_frags, round(total_frags - displacement)))
                        frame_ks.append((k, direction))
                    # reverse ramp: total_frags → 0
                    for i in range(rev_frames):
                        k = int(total_frags * (rev_frames - 1 - i) / (rev_frames - 1)) if rev_frames > 1 else 0
                        frame_ks.append((k, -1))
                else:
                    fwd_frames = args.frames
                    frame_ks = [
                        (int(total_frags * i / (fwd_frames - 1)) if fwd_frames > 1 else total_frags, 1)
                        for i in range(fwd_frames)
                    ]
                total_output_frames = len(frame_ks)

                os.makedirs(args.output_dir, exist_ok=True)
                base      = os.path.basename(path)
                res_tag   = f"_{tw}x{th}"
                top_tag   = f"_top{top_n}" if top_n is not None else ""
                fade_tag  = "_fade" if args.animate_fade else ""
                inv_tag   = "_invert" if args.animate_invert else ""
                out_path  = os.path.join(
                    args.output_dir,
                    f"{base}_circular{res_tag}{top_tag}_reveal{fade_tag}{inv_tag}_{uuid.uuid4()}.mp4",
                )
                enc_w = tw + (tw % 2)
                enc_h = th + (th % 2)
                cmd = [
                    "ffmpeg", "-y",
                    "-f", "image2pipe", "-vcodec", "png", "-r", str(args.fps), "-i", "-",
                    "-vcodec", "libx264", "-pix_fmt", "yuv420p",
                    "-vf", f"scale={enc_w}:{enc_h}",
                    out_path,
                ]
                proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)
                for frame_idx, (k, direction) in enumerate(frame_ks):
                    frame_frags = sorted_frags[:k]
                    print(f"Rendering frame {frame_idx + 1}/{total_output_frames} ({k}/{total_frags} fragments)...")
                    if args.animate_fade and k > 0:
                        # Gaussian falloff centred on the "active" end of the
                        # fragment list.  During forward pass (direction=1) the
                        # newest fragment is at index k-1; during the reverse
                        # pass (direction=-1) the fragment about to leave is at
                        # index 0, so we mirror the distance calculation.
                        idx_arr  = np.arange(k, dtype=np.float64)
                        dist     = (k - 1) - idx_arr
                        t        = dist / taper_window_count
                        bvals    = np.maximum(
                            np.exp(-args.taper_sharpness * t * t), 0.05
                        ).astype(np.float32)
                    else:
                        bvals = None
                    surface = render(frame_frags, n, progress=True, circle_fill=render_fill,
                                     brightness=bvals, draw_arcs=not args.no_arcs)
                    if (tw, th) != (_RENDER_SIZE, _RENDER_SIZE):
                        surface = downscale(surface, tw, th)
                    buf = io.BytesIO()
                    surface.write_to_png(buf)
                    proc.stdin.write(buf.getvalue())
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    print(f"ffmpeg exited with code {proc.returncode}", file=sys.stderr)
                    sys.exit(1)
                print(f"Saved: {out_path}")
                webbrowser.open(f"file://{os.path.abspath(out_path)}")
            else:
                do_circular = args.circular or not args.linear
                do_linear   = args.linear
                os.makedirs(args.output_dir, exist_ok=True)
                base = os.path.basename(path)
                res_tag = f"_{tw}x{th}"
                top_tag = f"_top{top_n}" if top_n is not None else ""
                arc_tag = "_noarcs" if args.no_arcs else ""
                if do_circular:
                    surface = render(render_frags, n, progress=True, circle_fill=render_fill, draw_arcs=not args.no_arcs)
                    if (tw, th) != (_RENDER_SIZE, _RENDER_SIZE):
                        surface = downscale(surface, tw, th)
                    out_path = os.path.join(args.output_dir, f"{base}_circular{res_tag}{top_tag}{arc_tag}_{uuid.uuid4()}.png")
                    surface.write_to_png(out_path)
                    print(f"Saved: {out_path}")
                    if args.jpg:
                        _save_jpg(surface, out_path)
                    webbrowser.open(f"file://{os.path.abspath(out_path)}")
                if do_linear:
                    surface = render_linear(render_frags, n, progress=True, tw=tw, th=th, draw_arcs=not args.no_arcs)
                    if (tw, th) != (_RENDER_SIZE, _RENDER_SIZE):
                        surface = downscale(surface, tw, th)
                    out_path = os.path.join(args.output_dir, f"{base}_linear{res_tag}{top_tag}{arc_tag}_{uuid.uuid4()}.png")
                    surface.write_to_png(out_path)
                    print(f"Saved: {out_path}")
                    if args.jpg:
                        _save_jpg(surface, out_path)
                    webbrowser.open(f"file://{os.path.abspath(out_path)}")
