BEGIN {
  OFS = "\t";
} {
  if (len) {
    if ($0~/[ ]*length=/) {
      split($0, a, "[ ]*length=");
      print a[1], a[2];
    }
  } else {
    if ($0~/[ ]*tmscore=/) {
      split($0, a, "[ ]*tmscore=");
      print a[1], a[2];
    }
  }
} END {
}
