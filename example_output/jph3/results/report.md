# JPH3: junction sequencing summary

jseq 1.1.0. Config `config/jph3.json`.

## JPH3_5p (spacing mode, tile filter: none)

- flowcell/lane `000H5VG57:1`: kept 24 of 24 tiles

Tile QC: all libraries PASS

## JPH3_3p (spacing mode, tile filter: none)

- flowcell/lane `000H5VG57:1`: kept 24 of 24 tiles

Tile QC: all libraries PASS

## Net-zero by culture
Culture rows: mean of that culture's libraries, SD between them. ALL row: mean, SD and SEM of the culture means, so each culture counts once however many libraries it has.

| amplicon | culture | libraries | mean % | SD | SEM |
|---|---|---|---|---|---|
| JPH3_5p | 1 | 1 | 70.21 |  |  |
| JPH3_5p | 2 | 1 | 74.78 |  |  |
| JPH3_5p | 3 | 1 | 75.12 |  |  |
| JPH3_5p | ALL | 3 cultures | 73.37 | 2.7426 | 1.5834 |
| JPH3_3p | 1 | 1 | 66.97 |  |  |
| JPH3_3p | 2 | 1 | 68.03 |  |  |
| JPH3_3p | 3 | 1 | 69.14 |  |  |
| JPH3_3p | ALL | 3 cultures | 68.04 | 1.0838 | 0.6257 |

Heatmap: `example_output/jph3/results/heatmap_by_culture.png` / `.pdf`; values in `example_output/jph3/results/heatmap_by_culture.csv`.
