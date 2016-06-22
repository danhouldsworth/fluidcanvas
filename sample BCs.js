function objectBoundary(ux, uy){

    var s4th        = Math.floor(SIZE / 4),
        s10th       = Math.floor(SIZE / 10);
    // -- SQUARE
    for (var y = s4th - s10th ; y < s4th + s10th; y++) {
        for (var x = s4th - s10th; x < s4th + s10th; x++) {
            ux[arrayIndex(x, y)] = 0;
            uy[arrayIndex(x, y)] = 0;
        }
    }
    // -- CIRCLE
    for (var y = SIZE - s4th - s10th; y < SIZE -s4th + s10th; y++) {
        for (var x = s4th - s10th; x < s4th + s10th; x++) {
            if ((x-s4th)*(x-s4th)+(y - SIZE + s4th)*(y - SIZE + s4th) < (s10th*s10th) ) {
                ux[arrayIndex(x, y)] = 0;
                uy[arrayIndex(x, y)] = 0;
            }
        }
    }

    // -- FLAT PLATE @ 30deg AOA
    // for (var y = SIZE/10; y < SIZE/5; y++) {
    //     ux[arrayIndex(y, SIZE *3.5/ 10 +0)] = 0;
    //     uy[arrayIndex(y, SIZE *3.5/ 10 +0)] = 0;
    //     ux[arrayIndex(y, SIZE *3.5/ 10 +1)] = 0;
    //     uy[arrayIndex(y, SIZE *3.5/ 10 +1)] = 0;
    //     ux[arrayIndex(y, SIZE *3.5/ 10 +2)] = 0;
    //     uy[arrayIndex(y, SIZE *3.5/ 10 +2)] = 0;
    //     ux[arrayIndex(y, SIZE *3.5/ 10 +3)] = 0;
    //     uy[arrayIndex(y, SIZE *3.5/ 10 +3)] = 0;
    //     ux[arrayIndex(y, SIZE *3.5/ 10 +4)] = 0;
    //     uy[arrayIndex(y, SIZE *3.5/ 10 +4)] = 0;
    // }
    // for (var y = SIZE/10; y < SIZE/5; y++) {
    //     ux[arrayIndex(2*y+0, SIZE/ 4 + y)] = 0;
    //     uy[arrayIndex(2*y+0, SIZE/ 4 + y)] = 0;
    //     ux[arrayIndex(2*y+1, SIZE/ 4 + y)] = 0;
    //     uy[arrayIndex(2*y+1, SIZE/ 4 + y)] = 0;
    //     ux[arrayIndex(2*y+2, SIZE/ 4 + y)] = 0;
    //     uy[arrayIndex(2*y+2, SIZE/ 4 + y)] = 0;
    //     ux[arrayIndex(2*y+3, SIZE/ 4 + y)] = 0;
    //     uy[arrayIndex(2*y+3, SIZE/ 4 + y)] = 0;
    //     ux[arrayIndex(2*y+4, SIZE/ 4 + y)] = 0;
    //     uy[arrayIndex(2*y+4, SIZE/ 4 + y)] = 0;
    // }
    // for (var y = SIZE/10; y < SIZE/5; y++) {
    //     ux[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-0)] = 0;
    //     uy[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-0)] = 0;
    //     ux[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-1)] = 0;
    //     uy[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-1)] = 0;
    //     ux[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-2)] = 0;
    //     uy[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-2)] = 0;
    //     ux[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-3)] = 0;
    //     uy[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-3)] = 0;
    //     ux[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-4)] = 0;
    //     uy[arrayIndex(SIZE * 3 / 10 + y, SIZE / 4 + 2*y-4)] = 0;
    // }

    // -- SHARP APERTURE
    // for (var y = 0; y < SIZE / 2 - 5; y++) {
    //     ux[arrayIndex(1 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(1 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(2 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(2 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(3 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(3 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(4 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(4 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(5 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(5 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(6 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(6 + SIZE / 4, y)] = 0;
    // }
    // for (var y = SIZE/2 + 5; y < SIZE; y++) {
    //     ux[arrayIndex(1 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(1 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(2 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(2 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(3 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(3 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(4 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(4 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(5 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(5 + SIZE / 4, y)] = 0;
    //     ux[arrayIndex(6 + SIZE / 4, y)] = 0;
    //     uy[arrayIndex(6 + SIZE / 4, y)] = 0;
    // }

}
