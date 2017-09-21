create table fits_and_extracted_data (
    x    INTEGER NOT NULL,
    y    INTEGER NOT NULL,
    fit_id INTEGER NOT NULL,
    successful_fit INTEGER NOT NULL,
    fit_attempt INTEGER NOT NULL,
    reduc_chisq REAL,
    pf   TEXT,
    pferr TEXT,
    p0 TEXT,
    fit_bounds TEXT,
    fwhm_data TEXT,
    PRIMARY KEY( x, y, fit_id )
);


