
## [1.1.3] - (16/04/2026)

- Setup of auto-versioning
## [1.1.2] - (14/04/2026)

- Fixing example 08/time-slices
- Fixing the warning due to the use of the deprecated datetime.utcnow() in the default filenames
- The default file format for pictures is now svg instead of png in save_time_slices()

## [1.1.1] - (25/12/2025)

- Adding AUTHORS, LICENCE files
- Adding "Citation" section in both documentation and README

## [1.1.0] - (04/12/2025)

- Paper.md now detail both the models and numerical schemes.
- The README.md and documentation usage.rst now present how to create a simple domain, mesh and how to run simple simulations.
- python version (>=3.10) --> (>=3.11) because of the use of **typing.Self** for type hint of special methods *__add__* and *__mult__* of CellValueMap and VertexValueMap objects.

## [1.0.0] - (27/10/2025)

_First release._
