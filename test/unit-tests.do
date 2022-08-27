version 14.1

* qui do test/unit-tests-consistency.do
qui do test/unit-tests-basic.do

capture program drop main
program main

    matrix l_vec0 = 0 \ 0 \ 0 \ 0
    matrix l_vec1 = 1 \ 1 \ 1 \ 1
    matrix l_vec2 = 0 \ 3 \ -2 \ 0
    matrix l_vec3 = 0 \ 0 \ 0 \ 0 \ 0 \ 0 \ 1
    matrix l_vec4 = -1 \ 2 \ -3 \ 4 \ -5 \ 6 \ -7

    unit_test basic_failures

    unit_test basic_checks, norelmag reference(4) method(FLCI)
    unit_test basic_checks, norelmag reference(4) method(C-F)
    forvalues i = 0 / 2 {
        unit_test basic_checks, norelmag reference(4) l_vec(l_vec`i') method(FLCI)
        unit_test basic_checks, norelmag reference(4) l_vec(l_vec`i') method(C-F)
    }
    unit_test basic_checks, norelmag reference(1) method(FLCI)
    unit_test basic_checks, norelmag reference(1) method(C-F)
    unit_test basic_checks, norelmag reference(7) method(FLCI)
    unit_test basic_checks, norelmag reference(7) method(C-F)
    forvalues i = 3 / 4 {
        unit_test basic_checks, norelmag reference(1) l_vec(l_vec`i') method(FLCI)
        unit_test basic_checks, norelmag reference(1) l_vec(l_vec`i') method(C-F)
    }

    foreach rm in norelmag relmag {
        unit_test basic_checks, `rm' reference(4) method(Conditional)
        unit_test basic_checks, `rm' reference(4) method(C-LF)
        forvalues i = 0 / 2 {
            unit_test basic_checks, `rm' reference(4) l_vec(l_vec`i') method(Conditional)
            unit_test basic_checks, `rm' reference(4) l_vec(l_vec`i') method(C-LF)
        }
        unit_test basic_checks, `rm' reference(1) method(Conditional)
        unit_test basic_checks, `rm' reference(1) method(C-LF)
        unit_test basic_checks, `rm' reference(7) method(Conditional)
        unit_test basic_checks, `rm' reference(7) method(C-LF)
        forvalues i = 3 / 4 {
            unit_test basic_checks, `rm' reference(1) l_vec(l_vec`i') method(Conditional)
            unit_test basic_checks, `rm' reference(1) l_vec(l_vec`i') method(C-LF)
        }
    }

    * NB: checks l_vec1 and l_vec2 took forever and a day
    matrix l_vec0 = J(23, 1, 0)
    matrix l_vec1 = J(23, 1, 1)
    mata st_matrix("l_vec2", ((mod(1::23, 2) :* 2) :- 1) :* (1::23))

    unit_test reg_checks, norelmag method(FLCI)
    unit_test reg_checks, norelmag method(C-F)
    forvalues i = 0 / 2 {
        unit_test reg_checks, norelmag l_vec(l_vec`i') method(FLCI)
        unit_test reg_checks, norelmag l_vec(l_vec`i') method(C-F)
    }
    foreach rm in norelmag relmag {
        unit_test reg_checks, `rm' method(Conditional)
        unit_test reg_checks, `rm' method(C-LF)
        forvalues i = 0 / 2 {
            unit_test reg_checks, `rm' l_vec(l_vec`i') method(Conditional)
            unit_test reg_checks, `rm' l_vec(l_vec`i') method(C-LF)
        }
    }
end

capture program drop unit_test
program unit_test
    syntax namelist(min=1 max=1), [NOIsily tab(int 4) *]
    cap `noisily' `namelist', `options'
    if ( _rc ) {
        di as error _col(`=`tab'+1') `"test(failed): `namelist', `options'"'
        exit _rc
    }
    else di as txt _col(`=`tab'+1') `"test(passed): `namelist', `options'"'
end

main
