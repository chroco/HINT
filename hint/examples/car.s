log "CASE car"
set rule list

buying in {1,2,3,4}
maint in {1,2,3,4}
doors in {1,2,3,4}
persons in {1,2,3}
lug_boot in {1,2,3}
safety in {1,2,3}

y in {1,2,3,4}

y depends on {buying, maint, doors, persons, lug_boot, safety}

car depends on {price, tech}
price depends on {buying, maint}
tech depends on {safety, comfort}
comfort depends on {doors, persons, lug_boot}

car,price,tech,comfort in {1,2,3,4}

sel y

load rule car_data using 1
write dot car to aaa

load act.s
