#!/usr/bin/awk -f
### Script to extract average pairwise identities from T-COFFEE output
### example to run: ./APIextract.awk <path to T-COFFEE Indentity file> > <desired output file>

BEGIN{
  ### print column description
  printf "Promoter,Indentity,NumberOfIsolates\n";
  N=0;
}
{
  ### extract and print promoter name
  if (($0 ~ /^#/) && ($3 ~ /^U00096.3_MG1655_genome/)) {
    split($3, P, "_");
    N=1;
  } else if ($0 ~ /^#/) {
    N+=1;
  }
  ### print total percent identity (API) of the sequence among strains
  if (($1 ~ /^TOT/)) {
    printf "%s,%s,%s\n", P[4], $4, N-1;
  }

}
