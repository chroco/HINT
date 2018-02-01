log "CASE nursery"

y depends on {parents,has_nurs,form,childs,housing,finance,social,health}

parents in {1,2,3}
has_nurs in {1,2,3,4,5}
form in {1,2,3,4}
childs in {1,2,3,4}
housing in {1,2,3}
finance in {1,2}
social in {1,2,3}
health in {1,2,3}

y in {1,2,3,4,5}

nursery in {1,2,3,4,5}
employ in {1,2,3,4}
struct_finan, structur, soc_health in {1,2,3}

nursery depends on {employ,struct_finan,soc_health}
employ depends on {parents,has_nurs}
struct_finan depends on {housing,finance,structur}
structur depends on {form,childs}
soc_health depends on {social,health}

#write dot nursery to aaa
#quit

sel y
load rule nursery.data using 1
load act.s
quit
