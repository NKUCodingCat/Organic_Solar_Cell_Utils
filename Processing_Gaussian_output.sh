cat v.out* > tmp.out
bash grp

cat E.out-*-Typ_2 | sed "s/\([0-9]\+\)m.\+\.log/\1/g" | awk '!a[$1]++' | awk '{printf "%s\t%s\t%s\r\n" , $1, $3, $5}' > Typ2_HOMO_LUMO.txt
cat E.out-*-Typ_1 | sed "s/\([0-9]\+\)m.\+\.log/\1/g" | awk '!a[$1]++' | awk '{printf "%s\t%s\t%s\r\n" , $1, $3, $5}' > Typ1_HOMO_LUMO.txt
#   cat all outputs | replace file name              | delete duplicate |     grep energy