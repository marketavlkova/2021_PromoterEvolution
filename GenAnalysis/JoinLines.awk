#!/usr/bin/awk -f
### Script for joining the lines of identical sequences from .aln alignment files
### example to run: ./JoinLines.awk <input file> > <output file>

{
  ID=$1;
  $1="";
  SEQ[ID]=SEQ[ID]$0;
}
END{
  printf "?\n";
  for (ID in SEQ) {
    printf "%s:%s\n", ID, SEQ[ID];
  }
  printf "!\n";
}
