# PrePH
Predict PanHandles
PrePH is a set of Python scripts for finding Pairs of Complementary regions 

Forked from [kalmSveta/PrePH](https://github.com/kalmSveta/PrePH) by Svetlana Kalmykova.

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Usage

This software predicts RNA structures between conserved regions within genes. It takes the following genome annotation as input:

+ `fasta` sequence;
+ `gtf` annotation (gzip compressed);
+ UCSC-formatted table with conserved elements (gzip compressed);
+ `bed` files with repeats to exclude (gzip compressed);

and returns a list of structures in tabular and `bigBed` format.

The software is implemented as Snakemake pipeline. To run it, set the parameters in the `config/config.yaml` file and run the following command:

```
snakemake --use-conda --configfile config/config.yaml -c8 all
```

Please refer to [PrePH README](workflow/scripts/PrePH/README.md) for the documentation for individual PrePH commands.