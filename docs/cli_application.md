# Roofer CLI App

{program}`roofer` is a command line application for generating LoD 1.2, LoD 1.3 and LoD 2.2 building models from pointclouds in the [CityJSONSeq](https://www.cityjson.org/cityjsonseq/) format.

With the {program}`roofer` CLI application you can:

+ Generate 3D building models for large areas with a reasonably small memory footprint. This is achieved by processing the input data in square tiles. No pre-tiling is needed.
+ Benefit from multi-threading for faster processing.
+ Read common vector input sources and filter them based on attributes or spatial extent. This means you do not need to crop the input data beforehand.
+ Obtain detailed information per building object about the reconstruction process and the quality of the generated models and the input point cloud.


## Usage
```{code-block} text
Automatic LoD 2.2 building reconstruction from airborne lidar pointclouds

Usage:
  roofer [options] <pointcloud-path>... <polygon-source> <output-directory>
  roofer [options] (-c | --config) <config-file> [(<pointcloud-path>... <polygon-source>)] <output-directory>
  roofer -h | --help
  roofer -H | --help-all
  roofer -a | --attributes
  roofer -v | --version
```

### Examples

```{code-block} shell
roofer pointcloud.laz footprints.gpkg output-dir
roofer --lod12 --lod22 pointcloud.laz footprints.gpkg output-dir
roofer --filter 'identificatie=1980100000265200' pointcloud.laz footprints.gpkg output-dir
roofer -c config.toml output-dir
roofer --tiling pointcloud.laz footprints.gpkg - | cjseq collect > buildings.city.json
```

### Positional arguments

```{option} <pointcloud-path>
Path to pointcloud file (.LAS or .LAZ) or folder that contains pointcloud files.
```

```{option} <polygon-source>
Path to roofprint polygon source. Can be an OGR supported file (eg. GPKG) or database connection string.
```

```{option} <output-directory>
Supply a directory for file output, or `-` to stream CityJSONSeq to stdout. The same value can be set with `output-directory = "-"` under `[output]` in TOML.
```

```{include} cli-options.md
```

## Output format
Roofer writes [CityJSONSequence](https://www.cityjson.org/cityjsonseq/) as newline-delimited JSON. Each feature represents one building, with an optional per-tile `TINRelief` terrain feature when `--terrain` is enabled.

> [!WARNING]
> Terrain output is experimental. Its geometry, attributes, file naming, and configuration may change, and the feature may be removed in a future release.

For terrain output, `--terrain-nodata-mode` controls whether only complete grid quads are triangulated, three-sample local triangles are retained, or isolated missing samples are interpolated before local triangulation.

> [!TIP]
>  In the near future we expect to introduce a converter program to convert the CityJSONSequence files to other formats like GPKG, OBJ and others.

### Streaming to stdout

With `-` as the destination, Roofer writes one CityJSON header and flushes each
building feature as its reconstruction finishes. A reader can process the first
building while its neighbors are still running. Console diagnostics go to stderr.
For example, this validates the features as they arrive.

```shell
roofer -c config.toml - | cjval
```

With `--tiling`, buildings from different tiles can arrive in any order. Every
tile uses the scale and translation in the single stream header. By default the
translation is the center of the processing region in XY and zero in Z. Use
`--cj-translate` to supply an explicit origin. An empty selection produces just the
header.


### Output attributes
```{include} output-attributes.md
```

### Underpass attributes

An explicit `underpass-height-attribute` value is used as an absolute ceiling
elevation. If the attribute is omitted from the configuration or its value is
null or unset, Roofer uses the building floor elevation plus
`underpass-default-height` (4 m by default). Only these fallback ceilings can
be lowered to fit below the roof. Lowering uses `underpass-roof-margin` and
stops at `underpass-min-height` above the highest ground in the passage. Where
the roof is at or below the resulting ceiling, the opening cuts through the
roof.

The building attribute `rf_has_underpass` is a boolean indicating whether the
reconstructed building contains at least one underpass ceiling. It is false
when all matching inputs are skipped or lie outside the reconstructed footprint,
and for LoD 1 and fallback extrusions. Rename or omit it using
`[output.attributes] has_underpass`, like the other output attributes.

Select attributes from the underpass polygon source to copy onto the matching
`OuterCeilingSurface` semantic objects:

```toml
[input]
underpass-source = "underpasses.gpkg"
underpass-attributes = ["underpass_id", "name"]
```

The CLI equivalent is `--underpass-attributes underpass_id,name` (a single name
is also accepted). Every selected name must exist in the source layer; missing
fields cause an error before reconstruction. Boolean, integer, real and string
values retain their types, and null or unset values remain null. Other OGR
field types use their text representation. Original attribute names are kept;
the reserved CityJSON `type` field remains `OuterCeilingSurface`.

Where polygons overlap, the input with the highest adjusted ceiling supplies
the attributes. Equal adjusted heights use the first input in source order.
Selecting attributes retains ceiling boundaries between different inputs,
including ceilings at the same height.

### Roof lines

`line-source` adds 2D roof lines, for example ridges and dormer outlines
extracted from orthophotos, to the lines roofer detects in the point cloud:

```toml
[input]
line-source = "roof-lines.gpkg"

[reconstruction.line-regulariser]
external-line-priority = 0
```

The source holds LineString or MultiLineString features in the footprint CRS;
other geometry types are skipped. For each building, the parts of the lines
inside its footprint are regularised together with the detected boundary and
plane intersection lines and then partition the roof. The graph cut decides
which of the resulting faces remain. `external-line-priority` sets their weight
in the regularisation, where each cluster of nearby parallel lines is replaced
by its highest priority lines:

- 0 (default): as detected boundary lines, averaged with those nearby and
  replaced by a nearby plane intersection line.
- 1: replace nearby boundary lines, replaced by a nearby plane intersection
  line.
- 2: as plane intersection lines. Unless `merge-intersection-lines` is set,
  these are kept separately, so a line close to a detected ridge adds a narrow
  face that the graph cut assigns to either side.

`roofprint-skeleton` adds the footprint's own candidate lines without a line
source: the arcs of its interior straight skeleton of at least 1 m, which are
the ridges and hips of a hipped roof over the footprint. They join the external
lines with the same priority. This helps where the point cloud is too smooth to
show ridges, for example a DSM from image matching, and costs about 15 % more
reconstruction time.

```toml
[reconstruction.line-regulariser]
roofprint-skeleton = true
```

### Narrow structures

A firewall or parapet on a flat roof is often narrower than twice the alpha
shape's radius (r = sqrt(`alpha-shaper.alpha`)). The region grower leaves its
points unsegmented or fits them as a narrow plane whose alpha shape breaks up,
so the roof face runs on underneath. `narrow-planes` adds a step after the
plane detector:

```toml
[reconstruction.plane-detector]
narrow-planes = true
narrow-max-slope = 10.0   # degrees, steepest roof a structure stands on
narrow-min-height = 0.5   # m over the roof on both sides
narrow-max-height = 3.0   # m over the highest roof plane's median
narrow-min-length = 4.0   # m, for structures away from the footprint
```

```{mermaid}
flowchart LR
  PD[plane detector] --> NP[narrow structures]
  NP -- "planes, labels" --> AS[alpha shaper]
  NP -- "rectangles" --> HF[heightfield]
  AS --> HF
  AS --> LR[line regulariser]
  NP -- "rectangle sides" --> ARR[arrangement]
  LR --> ARR --> OPT[graph cut]
  NP -- "rectangles" --> OPT
```

1. Points that stand more than `plane-epsilon` over a less steep plane within
   2r, unsegmented or in a plane that mostly stands, are split into straight
   runs at most sqrt(3)·r wide. A run that stands `narrow-min-height` over the
   roof on both sides, but not more than `narrow-max-height` over the
   median height of the building's highest plane (a wire, or the ridge of a roof the detector
   missed), becomes a plane, level across, outlined by its
   rectangle. Pieces in line at one height are merged; ends run on over
   unsegmented points, to the footprint within 4r, and across a structure they
   meet at a T or a corner.
2. Left out are structures along the footprint, shorter than 4r, and those
   shorter than `narrow-min-length` that end away from the footprint between
   roof of one height (chimneys, rooftop equipment). Their points go back to
   their plane; a building without a structure left is reconstructed as
   without the option.
3. Each rectangle goes into the heightfield; its plane's alpha shape gives no
   lines. Its sides, carried on by the regulariser's `extension` but no
   farther than it is long, go into the arrangement past the regulariser,
   which would merge them.
4. The graph cut allows a narrow plane's label only on faces that lie in its
   rectangle: elsewhere it costs more than any data and smoothness cost it
   could save. The faces beside it take the roof planes as usual.

Roof surfaces in a rectangle and on its plane carry the attribute `rf_narrow`.

## Config file
Below is an example of a [TOML](https://toml.io/en/) configuration file for the {program}`roofer` CLI application. It shows all the available options. Some of these options are also available as command line arguments, in case one option is provided both in the configuration file and as a command line argument, the command line argument takes precedence.

```{literalinclude} example_full.toml
:language: toml
```
