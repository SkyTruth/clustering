Tools and pipelines:

    timestamptosec INFILE OUTFILE

Converts timestamps in Y-M-D H:M:S UTC (string) format to seconds since epoch (integer)

    augment INFILE OUTFILE

Adds weight and sigma columns with default values (score and 1.0),
which allows the raw data to be used by the cluster animations for the
deepest zoom levels. OUTFILE is suitable as input to cluster
animation.

    timecluser INFILE OUTFILE

Clusters by time difference, within vessel id. OUTFILE is suitable as
input to cluster animation.

    split --outdir=DIRNAME INFILE

Splits a file by date into multiple files

    dbcluster INFILE OUTFILE

Clusters indata by distance using the DBSCAN algorithm. OUTFILE is
suitable as input to cluster animation.

    join INFILE1 INFILE2 ...INFILEN OUTFILE

Concatenates all INFILEX into OUTFILE. If INFILEX is suitable as input
to cluster animation then OUTFILE will be too.
