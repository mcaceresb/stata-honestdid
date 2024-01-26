* exit, clear
* stata14-mp

use debugging, clear
mata st_matrix("b", st_data(., "debugging1")')
mata st_matrix("V", st_data(., "debugging2-debugging14"))
matrix l_vec = 1 \ 0 \ 0 \ 0 \ 0
honestdid, pre(1/7) post(8/12) l_vec(l_vec) mvec(0.001) delta(sd) alpha(0.1) b(b) vcov(V)
