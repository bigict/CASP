BEGIN {
  FS = "/";
} {
  printf("%s\t%s\t%s/%s.pdb\n", $(NF-1), $0, d, $(NF-1));
} END {
}
