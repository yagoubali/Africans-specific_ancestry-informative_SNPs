#!/usr/bin/awk -f
#cut -f3,6,7,8,11  SGDP_metadata.279public.21signedLetter.44Fan.samples.txt >SGDP.panel
#sed -i 's/(//g' SGDP.panel
#sed -i 's/)//g' SGDP.panel
#cat  SGDP.panel | sed -e '1d' |  awk -f   SGDP_panel.awk - >SGDP_ALL.panel
BEGIN {OFS="\t"; print "sample","pop","super_pop","gender"} 
{
#if (NR==1) {next;}
 
if ($3=="Africa") {super_pop="AFR"}
else if ($3=="America") {super_pop="AMR"}
else if ($3=="SouthAsia") {super_pop="SAS"}
else if ($3=="EastAsia") {super_pop="EAS"}
else if ($3=="Oceania") {super_pop="OCE"}
else if ($3=="CentralAsiaSiberia" && $4=="China") {super_pop="EAS"}
else if ($3=="CentralAsiaSiberia" && $4=="Kyrgyzystan") {super_pop="WAS"}
else if ($3=="CentralAsiaSiberia"&& $4=="Russia") {super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Abkhazia"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Albania"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Armenia"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Bulgaria"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Czechosloviapre1989"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="England"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Estonia"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Finland"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="France"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Georgia"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Greece"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Hungary"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Iceland"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Iran"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Iraq"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Israel"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="IsraelCarmel"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="IsraelCentral"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="IsraelNegev"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Italy"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="ItalyBergamo"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="ItalySardinia"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="ItalyTuscany"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Jordan"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Norway"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="OrkneyIslands"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Poland"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Russia"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="RussiaCaucasus"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Spain"){ super_pop="EUR"}
else if ($3=="WestEurasia"  && $4=="Tajikistan"){ super_pop="SAS"}
else if ($3=="WestEurasia"  && $4=="Turkey"){ super_pop="WAS"}
else if ($3=="WestEurasia"  && $4=="Yemen"){ super_pop="WAS"}

#Sex code ('1' = male, '2' = female, '0' = unknown)
if($5=="M") {sex="male"}
else if($5=="F") {sex="female"}
else if($5=="U") {sex="unknown"}
{print $1,$2,super_pop,sex}

}

END {}




