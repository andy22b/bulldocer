#!/usr/bin/env bash
evname="$1"

cd ${evname}

ev_file=${evname}.txt
offset=0.2

p_file=${evname}_p_station_locs.txt
s_file=${evname}_s_station_locs.txt

p_ps=${evname}_p_station_map.ps
s_ps=${evname}_s_station_map.ps

evla=$(awk 'NR==2{print $8}' ${ev_file})
evlo=$(awk 'NR==2{print $9}' ${ev_file})

misc="-Rg -JA${evlo}/${evla}/90/10c"

gmt pscoast ${misc} -Dc -A1000 -Ggrey80 -W -P -K > ${p_ps}
echo "${evlo} ${evla} 0 360" | gmt psxy ${misc} -SW89.5d+a -Wthinner,black -P -K -O >> ${p_ps}
echo "${evlo} ${evla} 0 360" | gmt psxy ${misc} -SW80d+a -Wthinner,red -P -K -O >> ${p_ps}
echo "${evlo} ${evla} 0 360" | gmt psxy ${misc} -SW30d+a -Wthinner,red -P -K -O >> ${p_ps}

while read stlo stla az name
do
if (( $(echo "${az} < 180." | bc -l)  ))
then
text_angle=$(echo "90. - $az" | bc -l)
just=ML
else
text_angle=$(echo "270. - $az" | bc -l)
just=MR
fi

az_rad=$(echo "${az} * 3.141592654 /180" | bc -l)
x_offset=$(echo "${offset} * s (${az_rad})" | bc -l)
y_offset=$(echo "${offset} * c (${az_rad})" | bc -l)

cat << EOF | gmt psxy -W ${misc} -K -O >> ${p_ps}
${stlo} ${stla}
${evlo} ${evla}
EOF

echo "${stlo} ${stla} ${name}" | gmt pstext -F+jMC+f8 -D${x_offset}/${y_offset} ${misc} -K -O >> ${p_ps}



done < ${p_file}

awk '{print $1, $2}' ${p_file} | gmt psxy -St0.2 -Gdarkorange ${misc} -W -O >> ${p_ps}

gmt psconvert ${p_ps} -Tf -P -A
cd ..