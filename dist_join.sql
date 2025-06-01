CREATE TABLE points AS
    SELECT
        row_number() OVER () AS id,
        ST_POINT(x, y) as geom
    FROM
        range(0, 200, 2) AS r(x),
        range(0, 200, 2) AS r(y);

.timer on
SELECT * FROM points as l JOIN points as r ON ST_DWithin(l.geom, r.geom, 4);
