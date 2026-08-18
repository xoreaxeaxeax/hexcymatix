# hexcymatix

![errlook.exe](examples/errlook.exe_circular_3840x1080_top100000_de997a4e-5b62-467e-b120-bb6cda5eb6d1_invert_center.jpg)

`hexcymatix` is an exploratory research project, loosely envisioned as a
binary reverse engineering tool for extracting structural information from
repeated byte patterns in binary files.

The idea never materialized into a usable tool, but still yields some
fascinating representations of the striking self-symmetry in everyday data.

## Overview

Sometimes in reverse engineering you have an example of some data — a table, a
compressed blob, a protocol payload — and want to find similar data elsewhere in
a file, without yet being able to write a rule or parser that describes it.  You
want similarity without specification.

`hexcymatix` posits that repeated byte strings imply structure or relationship,
even if we don't know what that structure is.  By visualizing repetitions, we
can surface related data without knowing exactly what we're looking for.

For example, the Pensées passage below has several self-repeating sequences 
(' incapable of ', ' which he ', 'ing either ') — `hexcymatix` lays the file out
across two lines and connects any repeated byte sequences of at least 8 bytes:

![Pascal Pensees](examples/pascal_example.png)

The connections that emerge then hint at underlying structure and relationships
in the data. 

Repeating this process on an everyday binary file starts to automatically pull
out unexpected self-similarity, connections, and structure.  As an example, here
is the error lookup executable in a standard Windows install:

![errlook.exe](examples/errlook.exe_linear_3840x1080_top100000_noarcs_a46c130a-1e40-43ec-9b3d-95ee99cba520.jpg)

The most important part of the process is adding arcs between the fragment
connections to make the output look cool:

![errlook.exe](examples/errlook.exe_linear_3840x1080_top100000_ba6dccde-78a8-449b-a839-65b84c4aebee.jpg)

Having completely lost the plot, we can experiment with varied representations
of the data.  For example, the same file and connections rendered around the
perimeter of a circle appears as follows:

![errlook.exe](examples/errlook.exe_circular_3840x1080_top100000_de997a4e-5b62-467e-b120-bb6cda5eb6d1.jpg)

## Usage

```
python hexcymatix.py [options] file [file ...]
```

### Flags

| Flag | Description |
|------|-------------|
| `--circular` | Circular layout: byte offsets map to angles on a ring (default) |
| `--linear` | Parallel-bars layout: two horizontal bars connected by crossing quads |
| `--top N [N ...]` | Keep only the N longest fragments before rendering; pass several values to render one image per cutoff |
| `--no-arcs` | Draw chords only, skipping the arcs (substantially faster, less visual noise) |
| `--resolution RES` | Output resolution as `WxH` or a preset: `480p`, `720p`, `1080p`, `1440p`, `4k`, `uwfhd` (2560x1080), `suwfhd` (3840x1080) (default: `1080p`) |
| `--ultra` | Shorthand for `--resolution 4000x4000` |
| `--circle-fill F` | Fraction of the canvas diameter used by the circle (default: `0.25`, or `0.50` for the ultra-wide presets) |
| `--jpg` | Also write an 80% quality JPEG alongside each PNG |
| `--output-dir DIR` | Directory for output files (default: `output`) |

### Examples

```sh
# Circular visualization of an executable
python hexcymatix.py --circular guidgen.exe

# Linear visualization, top 500 fragments only
python hexcymatix.py --linear --top 500 guidgen.exe

# Batch, 4K output
python hexcymatix.py --circular --resolution 4k *.exe
```

## Output

### Executables

#### cmd.exe

![cmd.exe circular](examples/cmd.exe_linear_3840x1080_top15000_9fcb2dbb-e36f-4764-9550-ca92597ac5aa.jpg)

![cmd.exe circular](examples/cmd.exe_circular_3840x1080_top15000_bf9a30f2-4e58-4e10-82c1-1f2db644e28d.jpg)

#### calc.exe

![calc.exe linear](examples/calc.exe_linear_3840x1080_top15000_4a3edb69-2cb8-4491-9e7a-3562db0fd781.jpg)

![calc.exe circular](examples/calc.exe_circular_3840x1080_top15000_3ddaa20e-2d83-4d42-be1c-76b69348c84f.jpg)

#### explorer.exe

![explorer.exe linear](examples/explorer.exe_linear_3840x1080_top15000_fac84dbc-d1c4-48e0-8574-37f75b3918b1.jpg)

![explorer.exe circular](examples/explorer.exe_circular_3840x1080_top15000_a8a17a3f-b2ae-4e37-9de0-52f51c16a19d.jpg)

### Fonts

#### Arial TTF
![arial.ttf linear](examples/arial.ttf_linear_3840x1080_top15000_e3d2efc9-84a4-41f6-adea-a4b6fdce8221.jpg)

![arial.ttf circular](examples/arial.ttf_circular_3840x1080_top15000_099cd2bf-e1ac-4a41-bd8b-e7eb9a83f99a.jpg)

#### Comic Sans TTF
![comic.ttf linear](examples/comic.ttf_linear_3840x1080_top100000_b2e7a96e-abc6-46fb-9f88-7ba623f263bf.jpg)

![comic.ttf circular](examples/comic.ttf_circular_3840x1080_top100000_a703cbbe-e6c7-4575-a2e8-8564909ae87e.jpg)

#### Courier New TTF
![cour.ttf linear](examples/cour.ttf_linear_3840x1080_top15000_ab578f79-b87a-43a0-959a-c8cbae556dd7.jpg)

![cour.ttf circular](examples/cour.ttf_circular_3840x1080_top15000_b4e864fa-8cb9-433b-b50e-d1e0b9b0a5f2.jpg)

### Libraries

#### kernel32.dll

![kernel32.dll linear](examples/kernel32.dll_linear_3840x1080_top15000_d7434c90-0dfa-4ca5-86a1-32ec31d33c5f.jpg)

![kernel32.dll circular](examples/kernel32.dll_circular_3840x1080_top15000_cfffd064-b082-4f16-ac33-b4ee2f887eda.jpg)

#### spwizres.dll

![spwizres.dll linear](examples/spwizres.dll_linear_3840x1080_top100000_48342a63-3add-4f75-9ee7-ea5c0bc17bff.jpg)

![spwizres.dll circular](examples/spwizres.dll_circular_3840x1080_top100000_7c01b75c-f481-466d-b428-df03081ff30c.jpg)

#### user32.dll

![user32.dll linear](examples/user32.dll_linear_3840x1080_top15000_fc5e6a10-8d9a-4515-8269-7794660f1be7.jpg)

![user32.dll circular](examples/user32.dll_circular_3840x1080_top15000_29c3b21d-0cd6-4003-8f2b-76de10bb6794.jpg)

### Images

#### Windows 95 Default Background ("Clouds")

![Clouds.bmp linear](examples/Clouds.bmp_linear_3840x1080_top100000_96f2b4dd-6226-4594-b933-dfa433d153a6.jpg)

![Clouds.bmp circular](examples/Clouds.bmp_circular_3840x1080_top100000_2d8cb7ba-1dca-434b-bbc1-39b2bf9dfe30.jpg)

#### Windows XP Default Background ("Bliss")

![bliss.jpg linear](examples/bliss.jpg_linear_3840x1080_top100000_537429ae-a6c0-4119-a7ed-9eec57505e04.jpg)

![bliss.jpg circular](examples/bliss.jpg_circular_3840x1080_top100000_5e75909f-27e3-48dc-aa65-1d4557f72742.jpg)

### Plain Text

#### Windows XP EULA

![eula.txt linear](examples/eula.txt_linear_3840x1080_top15000_6472a6c3-a717-4b96-ab7b-d30034bdb963.jpg)

![eula.txt circular](examples/eula.txt_circular_3840x1080_top15000_13de0aa1-42bc-4bce-80e5-3a70f7a91ec4.jpg)

#### Pascal's Pensées

![pg18269.txt linear](examples/pg18269.txt_linear_3840x1080_top15000_89061902-a6a0-4b0e-9209-279c8a71d5dd.jpg)

![pg18269.txt circular](examples/pg18269.txt_circular_3840x1080_top15000_35a5424c-b6d9-447d-8d8b-b69c6a213fcc.jpg)

### Audio

#### Windows XP startup sound

![xpstartu.wav linear](examples/xpstartu.wav_linear_3840x1080_top100000_0e4acaf1-6c57-4be8-9f70-39dc32f0f024.jpg)

![xpstartu.wav circular](examples/xpstartu.wav_circular_3840x1080_top100000_b5afd6bd-6150-4b1c-9243-811e3db2a31c.jpg)

#### Windows XP WMA (Beethoven's 9th, included with install)

![beethov9.wma linear](examples/beethov9.wma_linear_3840x1080_top100000_f13b85f6-88b4-4f37-bab4-490da2980d56.jpg)

![beethov9.wma circular](examples/beethov9.wma_circular_3840x1080_top100000_cb9878a2-6b5f-4386-a698-cc1c49a7fa57.jpg)

## Dependencies

```
pip install pycairo numpy pillow
```

`pycairo` builds against the native cairo library, which has to be present
first — `brew install cairo pkg-config` on macOS, or `apt install
libcairo2-dev pkg-config` on Debian/Ubuntu.

## Applications and Future

The tool was developed as a reverse engineering aid: given a single known
example of some data — a struct, a compressed blob, a protocol payload — find
other instances of it in the file without needing additional RE.  The hope was
that regions sharing many repeated subsequences with the known example would
cluster into candidates worth investigating further.

Ultimately, the idea seemed shaky and the use-case narrow, so the technique was
never developed into the sort of seamless plugin it needed to be for useful
research.  The work is shared here so that it is not lost to the ether, in the
hopes that perhaps someone else can build on it.

## AI Disclaimer

This project was originally written by hand in C#.  It has been ported to Python
using LLMs, with little oversight or review.

## Author

`hexcymatix` is a research effort from Christopher Domas (@xoreaxeaxeax).

![errlook.exe](examples/errlook.exe_circular_3840x1080_top100000_de997a4e-5b62-467e-b120-bb6cda5eb6d1_invert.jpg)