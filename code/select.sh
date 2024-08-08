#/bin/sh -f
# example: 

ls output/*prob* > t1
cut -c 8-17 t1 > t2
grep -Fvf t2 tilelist > tilelist2
