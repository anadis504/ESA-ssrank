
#include <stdio.h>
#include <stdint.h>

//Should probably be simplified to reduce number of variables
void print_field(uint64_t x, uint64_t i){
   //uint32_t shift_unit = 8 + 8*(i/4); //8 or 16
   //uint32_t mask = (1<<shift_unit) - 1; //0b11111111 or 0b1111111111111111
   //uint32_t shifty = 3 + (i/4); //3 or 4
   //uint32_t local_idx = i%4; //4..5 -> 0..1; 0..3 -> 0..3
   //uint32_t shift = local_idx << shifty;
   uint32_t mask = ((1 << (8 + 8*(i/4))) - 1);
   uint32_t shift = (i%4) << (3 + (i/4));
   uint32_t v = ((uint32_t *)&x)[i/4];
   uint32_t s = ((v>>shift)&mask);
   fprintf(stderr,"s: %d\n",s);
}

int main(int argc, char **argv){
   uint64_t x;
   uint8_t *px8 = (uint8_t *)&x;
   uint16_t *px16 = (uint16_t *)&x;
   px8[0] = 5;
   px8[1] = 23;
   px8[2] = 56;
   px8[3] = 123;
   px16[2] = 567;
   px16[3] = 1111;
   
   for(int i=0;i<6;i++){
      print_field(x,(uint64_t)i);
   }
   
   return 0;
}
