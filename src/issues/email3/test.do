import delimited using "src/issues/email3/V.txt", clear

mata
v = st_data(., 1)
b = (         0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0),
    (         0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0),
    (         0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0),
    (         0,            0,            0,            0,            0,            0,   -.03004681,   -.38323671,   -.22811466,   -.37702492,   -.44935237,   -.44881702),
    (-.22187288,    .09401739,   -.42721282,            0,   -.72195146,   -.57850283,   -.61602494,   -.60325407,   -.56625665,   -.34880705,   -.85817923,   -.41421178),
    (-.61101216,   -.80646463,    .05407468,    .28656434,   -.07872118)
V = J(length(b), length(b), .)
k = 0
for(i = 1; i <= length(b); i++) {
    for(j = 1; j <= i; j++) {
        V[i, j] = V[j, i] = v[++k]
    }
}
st_matrix("eb", b)
st_matrix("eV", V)
end

The first entries are all 0s...they're also all omitted. He probably just needs to pass -omit- and it will work
honestdid, pre(1/9) post(11/20) mvec(0.5(0.5)2) b(eb) vcov(eV)
