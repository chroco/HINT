seed 8410
copy y to table tmp

split 30% of tmp to y and tmp1

log relieff
#log inform
log gini
log gain

set distribution
set color on

rm red inform y

set dec global
set dec 2 to 3
set dec cm
#set dec sdtic

dec y

sel y
ls struct
compare y to table tmp1

