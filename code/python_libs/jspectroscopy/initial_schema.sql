create table fits (
    x    INTEGER NOT NULL,
    y    INTEGER NOT NULL,
    fit_id INTEGER NOT NULL,
    successful_fit INTEGER NOT NULL,
    npeaks INTEGER NOT NULL,
    last_attempt INTEGER NOT NULL,
    params_guess BLOB,
    fit_bounds BLOB,
    peak_guesses BLOB,
    model BLOB,
    PRIMARY KEY( x, y, fit_id )
);



create table metadata ( 
    xdim INTEGER NOT NULL,
    ydim INTEGER NOT NULL,
    peak_types BLOB,
    yield_data BLOB,
    timestamp TEXT
);
