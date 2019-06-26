void intDigits(int * digits, int nInput, int base, int len) {
  int i, n=nInput;
  for(i=len-1; i>=0; i--){
    digits[i]= n % base;
    n/= base;
  }
}

