# EnergyCommunities

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

## Todo
- variable and fix costs
- CO2 emission
    - Probably, emission factors need to be commonly declare...
    - How can one deal with the allocation of emissions where there is an internal trade?


