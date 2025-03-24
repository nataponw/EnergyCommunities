CREATE TABLE id (
    agentid       TEXT,
    description   TEXT   DEFAULT "",
    location_lat  REAL   DEFAULT 0.0,
    location_long REAL   DEFAULT 0.0,
    PRIMARY KEY (
        agentid
    )
);
INSERT INTO id (agentid) VALUES ('dummypeer');

CREATE TABLE pSca (
    parameter TEXT,
    value     REAL DEFAULT -99.9,
    PRIMARY KEY (
        parameter
    )
);
INSERT INTO pSca VALUES ('thCur_baseprice', 0.0);
INSERT INTO pSca VALUES ('thCur_incrementprice', 0.0);
INSERT INTO pSca VALUES ('thCur_max', 0.0);
INSERT INTO pSca VALUES ('inttrade_fee_el', 0.0);
INSERT INTO pSca VALUES ('inttrade_fee_th', 0.0);
INSERT INTO pSca VALUES ('flexev_share', 0.0);
INSERT INTO pSca VALUES ('flexev_duration', 0.0);

CREATE TABLE pY (
    year        INTEGER,
    impprice_el REAL    DEFAULT 99.9,
    expprice_el REAL    DEFAULT 00.0,
    impprice_th REAL    DEFAULT 99.9,
    expprice_th REAL    DEFAULT 00.0,
    impprice_bm REAL    DEFAULT 99.9,
    impprice_ng REAL    DEFAULT 99.9,
    impprice_ho REAL    DEFAULT 99.9,
    impprice_h2 REAL    DEFAULT 99.9,
    PRIMARY KEY (
        year
    )
);

CREATE TABLE pTec (
    tec                   TEXT,
    pot                   REAL      DEFAULT 0.0,
    scap_cap              REAL      DEFAULT 0.0,
    eff_el                REAL      DEFAULT 0.0,
    eff_th                REAL      DEFAULT 0.0,
    eff_st                REAL      DEFAULT 0.0,
    eff_ch                REAL      DEFAULT 0.0,
    eff_aux               REAL      DEFAULT 0.0,
    cost_var              REAL      DEFAULT 0.0,
    cost_fix              REAL      DEFAULT 0.0,
    scap_finalyear        INTEGER   DEFAULT 9999,
    scap_payment          REAL      DEFAULT 0.0,
    scap_finalpaymentyear INTEGER   DEFAULT 9999,
    cost_dec               REAL      DEFAULT 0.0,
    ncap_loanduration     INTEGER   DEFAULT 0,
    ncap_loaninterest     REAL      DEFAULT 0.0001,
    ncap_lifetime         INTEGER   DEFAULT 0,
    PRIMARY KEY (
        tec
    )
);

CREATE TABLE pYTec (
    year INTEGER,
    tec  TEXT,
    inv  REAL    DEFAULT 0.0,
    PRIMARY KEY (
        year,
        tec
    )
);

CREATE TABLE pYTS (
    year            INTEGER,
    timestep        INTEGER,
    dem_el_hh       REAL    DEFAULT 0.0,
    dem_el_ev       REAL    DEFAULT 0.0,
    dem_th_hh       REAL    DEFAULT 0.0,
    gen_pv          REAL    DEFAULT 0.0,
    gen_wt          REAL    DEFAULT 0.0,
    cop             REAL    DEFAULT 1.0,
    PRIMARY KEY (
        year,
        timestep
    )
);
