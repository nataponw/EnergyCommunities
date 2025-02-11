# EnergyCommunities

This package presents an operation and investment optimization model for energy communities with the options for local energy (electricity and heat) trading.
The core model is documented in the paper [Link](www.google.com).

An example script together with databases is provided in the folder `example`. In the example, the status quo scenario (SQ) and the electricity trading scenario (EL), as defined in the original paper are calculated.

## Specific types
Custom types---`Year`, `Peer`, and `Tec`---are created to differentiate these sets from timesteps.
These types and their associated methods are defined in `src\SetStructs.jl`.

## Database
A database is dedicated for an individual participant. A blank database can be created from the file `src\auxiliary_src\db_structure.sql` via the function `createblankdatabase`. Each database contains the following tables:

- `id`, a table with an agent id and location information
- `pSca`, a table with scalar parameters, e.g., trading fees
- `pY`, " year-dependent parameters, e.g., energy prices
- `pTec`, " technology-dependent parameters, e.g., efficiency
- `pYTec`, " parameters dependent of years and technologies, e.g., investment costs
- `pYTS`, " parameters dependent of years and timesteps, e.g., demand- or generation profiles

## Acknowledgment

<img align="right" width="180" height="180" src="images/Logo_BMWK.png">

This work was conducted in the research project BKM 2.0, which was sponsored by the Federal Ministry for Economic Affairs and Climate Action of Germany (*Bundesministerium für Wirtschaft und Klimaschutz* BMWK) under the grant 03EI4015A.

## Licence

Copyright (C) 2025 Fraunhofer ISE
<img align="right" width="160" height="50" src="images/Logo_FhgISE.png">

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
