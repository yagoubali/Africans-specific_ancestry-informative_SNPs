#!/usr/bin/awk -f
BEGIN {OFS="\t"; print "sample","pop","super_pop","gender"} 
{
if (NR==1) {next;}
 
if ($7=="ASW") {super_pop="AFR"}
else if ($7=="CEU") {super_pop="EUR"}
else if ($7=="CHB") {super_pop="EAS"}
else if ($7=="CHD") {super_pop="EAS"}
else if ($7=="GIH") {super_pop="SAS"}
else if ($7=="JPT") {super_pop="EAS"}
else if ($7=="LWK") {super_pop="AFR"}
else if ($7=="MEX") {super_pop="AMR"}
else if ($7=="MKK") {super_pop="AFR"}
else if ($7=="TSI") {super_pop="EUR"}
else if ($7=="YRI") {super_pop="AFR"}

#Sex code ('1' = male, '2' = female, '0' = unknown)
if($5==1) {sex="male"}
else if($5==2) {sex="female"}
else if($5==0) {sex="unknown"}
{print $2,$7,super_pop,sex}

}

END {}
