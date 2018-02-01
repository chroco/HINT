seed 400
copy y to table tmp
split 50% of tmp to y and tmp1
30% noise to y
set apriory y

write y tmp1 to xxx

set distribution
set color on
#set dec global
set dec 2 
set dec m
set noise m error
set dec deb 2
set m 2

# cross validation for estimation of m
#cv copy y to table tmp2
#cv set cv 10 for tmp2
#cv split tmp2 to y and tmp3 using XXX
#cv set apriory y
#cv dec y
#cv compare y to table tmp3
#cv quit


dec y

ls struct

# statistics
compare y to table tmp1
compare struct y xxx
count leaf y
sdtic y
quit
