/*
simpleTunnel - Navier Stokes 2D solver on CPU / Canvas.
Aims :
1. Validate my understanding of this kind of FEA method - Done. I'm going for unbounded iterative approach on cartesian grid. They go for convergent solver on mesh grid.
2. Verify that this approach to NS gives real seeming solutions - YEP!! Vortex streets and pressure drag WoopWoop!!
3. Read the paper - and verify that pixelated cartesian grid seems fine vs complex custom mesh (that's what was holding me back from FEM before) - See 1 above. It's a lot simpler than they all made out.
4. Is square parcels okay v gaussian mask? - Yes, but now I understand iterative vs solver, I'm going to do one single pressure calc per iteration, and so a gaussian mask will be more appropriate.
5. What are the density, Mach implications? - See 4, we need to start mapping to Mach, density, viscosity. And so have them consistent in our calcs independent of grid size.
6. Is static pressure, dynamic pressure and temperature clear? - Not sure the physical concept of dynamic pressure is that helpful in the real world anyway. Temperature is linked to divergence and we need to think more on it, see 5 above.
7. What if boundary conditions had through put of flow? Wrapped / tunnel etc. - Ha! Have done this do do objects. And now I'm clear on iterative approach per modelParams.deltaT we don't need boundary conditions.
8. How to have object boundaries within - change the indexing from full raster, to exclusion of object and boundary? - Done. No slip, no thru-flow. Think about div / pressure calc. Think about compuationally effective method of skipping this area.
9. Simply coding to my style. Then optimise so still understand physics and parameters. - Done
10. Is viscosity and shear force clear? - Wasn't here!! Need to add carefully, thinking about sheer forces and point 5 above.
** How can we record / take measurements along the wing? - Easy, trace the path of object, taking pressure readings.

ADD function for returning shape boolean - Also think about how return surface, so can get pressure vector / friction vector into L / D / moment
FIX boundary conditions calc for div & pressure at boundary. Eg. Why to we get reflection back of pressure from tunnel exit? Why does pressure get through thin object?
ADD steakline view mode
ADD timeline view mode
Add vector arrow view mode

Add known shapes
ADD airfoil shape

Output - Lift / Drag / Moment of object

*/
"use strict";/* jshint
browser       : true
*/

var canvas  = document.getElementById('canvas'),
    ctx     = canvas.getContext('2d'),
    mouseSample = document.getElementById('mouseSample'),
    SIZE    = canvas.width = canvas.height = canvas.clientWidth,
    rafRef,
    imageData,
    replayFrames= [],
    replayIndex = 0,
    visualScale = [],
    visualOffset= [],
    visualisationMode = 1,
    realWorldParams = {                         // ** All in SI
        tunnelLength            : 5,            // m
        tunnelSpeed             : 15,           // m/s
        density                 : 1.292,        // kg/m3                        @ STP
        kinematicViscosity      : 0.0000133,    // m2/s                         @ STP
        dynamicViscosity        : 0.0000172,    // kg/m.s OR Pa.s OR N.s/m2     @ STP
        standardStaticPressure  : 101325,       // kilo Pa                      @ STP
        speedOfSound            : 330           // m/s                          @ STP
    },
    parcel = {                              // ** All in SI **
        size    : (realWorldParams.tunnelLength / SIZE),
        area    : (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE),
        volume  : (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE)
    },
    modelParams = {                             // *** Length / Area / Volume in Parcels *** Time in iteratinons ***
        speedOfSound: 1,
        tunnelSpeed : 1    * realWorldParams.tunnelSpeed / realWorldParams.speedOfSound,
        deltaT      : 1 / (SIZE * realWorldParams.speedOfSound / realWorldParams.tunnelLength) // (Real world fluid) seconds per iteration ~ 15micro seconds
    },
    gaussianMask = [
        [0.3,1.0,3.0],
        [1.0,0.0,1.0],
        [3.0,1.0,3.0]
    ],
    gaussianWt  = gaussianMask[0][0] + gaussianMask[1][0] + gaussianMask[2][0] + gaussianMask[0][1] + gaussianMask[1][1] + gaussianMask[2][1] + gaussianMask[0][2] + gaussianMask[1][2] + gaussianMask[2][2],
    mouseX = 0,
    mouseY = 0,
    trunc3dp = function(raw) {return (raw * 1000 | 0) / 1000;},
    sampleField = function(e) {
        mouseX = (e.clientX - canvas.getBoundingClientRect().left) | 0;
        mouseY = (e.clientY - canvas.getBoundingClientRect().top ) | 0;
        reportSample();
    },
    reportSample = function(){
        mouseSample.innerHTML = "<p>At mouse pointer :</p>";
        mouseSample.innerHTML += "<p>Pressure : " + trunc3dp(p0[arrayIndex(mouseX, mouseY)]) + " <strong>Pa</strong></p>";
        mouseSample.innerHTML += "<p>Flux : " + trunc3dp(div[arrayIndex(mouseX, mouseY)]) + " <strong>kg/s</strong> (" + trunc3dp(div[arrayIndex(mouseX, mouseY)]) + " px/calc)</p>";
        mouseSample.innerHTML += "<p>Vx : " + trunc3dp(u0x[arrayIndex(mouseX, mouseY)]) + " <strong>m/s</strong> (" + trunc3dp(u0x[arrayIndex(mouseX, mouseY)]) + " px/calc)</p>";
        mouseSample.innerHTML += "<p>Vy : " + trunc3dp(u0y[arrayIndex(mouseX, mouseY)]) + " <strong>m/s</strong> (" + trunc3dp(u0y[arrayIndex(mouseX, mouseY)]) + " px/calc)</p>";
    };
canvas.addEventListener('mousemove', sampleField);



gaussianMask[0][0] /= gaussianWt;
gaussianMask[1][0] /= gaussianWt;
gaussianMask[2][0] /= gaussianWt;
gaussianMask[0][1] /= gaussianWt;
gaussianMask[1][1] /= gaussianWt;
gaussianMask[2][1] /= gaussianWt;
gaussianMask[0][2] /= gaussianWt;
gaussianMask[1][2] /= gaussianWt;
gaussianMask[2][2] /= gaussianWt;

ctx.fillRect(0, 0, SIZE, SIZE);
imageData = ctx.getImageData(0, 0, SIZE, SIZE);

var arraySize       = SIZE * SIZE,
    arrayIndex      = function (x, y){
        if (x<0) x=0;
        if (x>SIZE-1) x=SIZE-1;
        if (y<0) y=0;
        if (y>SIZE-1) y=SIZE-1;
        return (x + y * SIZE);
    },
    p0  = new Float32Array(arraySize),
    p1  = new Float32Array(arraySize),
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    u0x = new Float32Array(arraySize),
    u0y = new Float32Array(arraySize),
    u1x = new Float32Array(arraySize),
    u1y = new Float32Array(arraySize),
    div = new Float32Array(arraySize);
    // *****


// Initial conditions
for(var i = 0; i < arraySize; i++) {
    p1[i]   = p0[i]     = realWorldParams.standardStaticPressure;
    u1x[i]  = u0x[i]    = 0;
    u1y[i]  = u0y[i]    = 0;
}
// --

function physics(){

    // State at t = t0
    tunnelBoundary(u0x, u0y);
    objectBoundary(u0x, u0y);
    // - these constrained fields will be used for render / and t1 calcs
    computeDivergence(); // ideally the BCs have been clever about setting div indirectly
    boylesIdealGas();
    // pressureDissapation();
    // Evolve state at t = t1 = t0 + modelParams.deltaT
    convectionOfFluid(u0x, u1x);
    convectionOfFluid(u0y, u1y);
    // Pressure forces and Shear forces are considered to have acted on parcel for whole of discretised time step
    // Hence we apply pressure / shear forces from t0 onto velocities at t0 to get velocities at t1
    // pressureGradientAcceleration();
    shearForceAcceleration();

    tunnelBoundary(u1x, u1y);
    objectBoundary(u1x, u1y);
    //

    // Double buffer back to
    var temp;
    temp = p0;
    p0  = p1;
    p1  = temp;
    temp = u0x;
    u0x = u1x;
    u1x = temp;
    temp = u0y;
    u0y = u1y;
    u1y = temp;
}

function tunnelBoundary(ux, uy) {
    for (var y = 0; y < SIZE; y++) {
        for (var x = 0; x < SIZE; x++) {
            ux[arrayIndex(x, y)] += 0.0001; // 500 iterations to get to flow speed
        }
    }

    for (var x = 0; x < SIZE - 1; x++) {
        ux[arrayIndex(x, 0)] = 0;
        uy[arrayIndex(x, 0)] = 0;
        ux[arrayIndex(x, 1)] = 0;
        uy[arrayIndex(x, 1)] = 0;
        ux[arrayIndex(x, 2)] = 0;
        uy[arrayIndex(x, 2)] = 0;
        ux[arrayIndex(x, SIZE - 1)] = 0;
        uy[arrayIndex(x, SIZE - 1)] = 0;
        ux[arrayIndex(x, SIZE - 2)] = 0;
        uy[arrayIndex(x, SIZE - 2)] = 0;
        ux[arrayIndex(x, SIZE - 3)] = 0;
        uy[arrayIndex(x, SIZE - 3)] = 0;
    }
    for (var y = 3; y < SIZE - 4; y++) {
        ux[arrayIndex(0, y)] = modelParams.tunnelSpeed;
        uy[arrayIndex(0, y)] = 0.0;
        ux[arrayIndex(1, y)] = modelParams.tunnelSpeed;
        uy[arrayIndex(1, y)] = 0.0;
        ux[arrayIndex(2, y)] = modelParams.tunnelSpeed;
        uy[arrayIndex(2, y)] = 0.0;
    }
}

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

function convectionOfFluid(src, dest) {

    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION

    function bilerp(x, y) {
        // x,y as fractional pixel index
        // Linear interpolation between 2 values a-->b, parametised by 0<c<1
        var lerp = function (a, b, c) {
            if (c<0 || c>1) throw c;
            return a * (1 - c) + b * c;
        };
        var x0 = Math.floor(x),
            y0 = Math.floor(y),
            xm = x - x0,
            ym = y - y0,
            x1 = x0 + 1,
            y1 = y0 + 1,
            p00 = src[arrayIndex(x0, y0)],
            p01 = src[arrayIndex(x0, y1)],
            p10 = src[arrayIndex(x1, y0)],
            p11 = src[arrayIndex(x1, y1)];
        return lerp(lerp(p00, p10, xm), lerp(p01, p11, xm), ym);
    }

    for (var y = 1; y < SIZE - 1; y++) {
        for (var x = 1, dx, dy; x < SIZE - 1; x++) {
            dx = u0x[arrayIndex(x, y)];                         // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
            dy = u0y[arrayIndex(x, y)];                         // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
            dest[arrayIndex(x, y)] = bilerp(x - dx, y - dy);    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
        }
    }
}

function computeDivergence() {
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    // Divergence = net flux of fluid (parcels / iteration)
    for (var y = 1; y < SIZE - 1; y++) {
        for (var x = 1; x < SIZE - 1; x++) {
            div[arrayIndex(x, y)]  = u0x[arrayIndex(x + 1, y)];
            div[arrayIndex(x, y)] -= u0x[arrayIndex(x - 1, y)];
            div[arrayIndex(x, y)] += u0y[arrayIndex(x,     y + 1)];
            div[arrayIndex(x, y)] -= u0y[arrayIndex(x,     y - 1)];
            // Assumes we don't need to track the different densities of where the fluid came from / going to.
        }
    }
}

function boylesIdealGas() {
    for (var y = 1; y < SIZE - 1; y++) {
        for (var x = 1; x < SIZE - 1; x++) {
            p1[arrayIndex(x, y)] = p0[arrayIndex(x, y)] * (1 - div[arrayIndex(x, y)]);
        }
    }
}

function pressureDissapation() {
    for (var y = 1; y < SIZE - 1; y++) {
        for (var x = 1; x < SIZE - 1; x++) {
            p1[arrayIndex(x, y)] += gaussianMask[0][0] * p0[arrayIndex(x-1, y-1)];
            p1[arrayIndex(x, y)] += gaussianMask[1][0] * p0[arrayIndex(x+0, y-1)];
            p1[arrayIndex(x, y)] += gaussianMask[2][0] * p0[arrayIndex(x+1, y-1)];
            p1[arrayIndex(x, y)] += gaussianMask[0][1] * p0[arrayIndex(x-1, y+0)];
            p1[arrayIndex(x, y)] += gaussianMask[1][1] * p0[arrayIndex(x+0, y+0)];
            p1[arrayIndex(x, y)] += gaussianMask[2][1] * p0[arrayIndex(x+1, y+0)];
            p1[arrayIndex(x, y)] += gaussianMask[0][2] * p0[arrayIndex(x-1, y+1)];
            p1[arrayIndex(x, y)] += gaussianMask[1][2] * p0[arrayIndex(x+0, y+1)];
            p1[arrayIndex(x, y)] += gaussianMask[2][2] * p0[arrayIndex(x+1, y+1)];
        }
    }
}

function pressureGradientAcceleration() {
    var mass        = realWorldParams.density * parcel.volume;
    var realToModel = modelParams.deltaT / parcel.size; // seconds per ITER / meters per parcel ==> real speed (m/s) to parcel/iteration
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    for (var y = 1; y < SIZE - 1; y++) {
        for(var x = 1, px0, px1, py0, py1, gradPx, gradPy, Fx, Fy; x < SIZE - 1; x++) {
            px0 = p0[arrayIndex(x - 1,   y)];
            px1 = p0[arrayIndex(x + 1,   y)];
            py0 = p0[arrayIndex(x,       y - 1)];
            py1 = p0[arrayIndex(x,       y + 1)];
            gradPx = (px1 - px0) / 2; // 2 cells wide from centre of sample (n-1) to centre of sample (n+1)
            gradPy = (py1 - py0) / 2;
            // Increasing pressure gradient deccelerates fluid
            Fx = -gradPx * realWorldParams.parcelArea ;
            Fy = -gradPy * realWorldParams.parcelArea ;
            u1x[arrayIndex(x, y)] += realToModel * Fx / mass;
            u1y[arrayIndex(x, y)] += realToModel * Fy / mass;
        }
    }
}

function shearForceAcceleration(){
    var mass        = realWorldParams.density * parcel.volume;
    var realToModel = modelParams.deltaT / parcel.size; // seconds per ITER / meters per parcel ==> real speed (m/s) to parcel/iteration
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    for (var y = 1; y < SIZE - 1; y++) {
        for(var x = 1, sheerX0, sheerX1, sheerX2, sheerY0, sheerY1, sheerY2, Fx, Fy; x < SIZE - 1; x++) {
            sheerY0 = u0y[arrayIndex(x - 1, y)];
            sheerY1 = u0y[arrayIndex(x,     y)];
            sheerY2 = u0y[arrayIndex(x + 1, y)];
            sheerX0 = u0x[arrayIndex(x,     y - 1)];
            sheerX1 = u0x[arrayIndex(x,     y)];
            sheerX2 = u0x[arrayIndex(x,     y + 1)];
            // dynamicViscosity applies Force on fluid parcel in direction of strain force
            Fx  = (sheerX0 - sheerX1) * realWorldParams.dynamicViscosity * parcel.area;
            Fx += (sheerX2 - sheerX1) * realWorldParams.dynamicViscosity * parcel.area;
            Fy  = (sheerY0 - sheerY1) * realWorldParams.dynamicViscosity * parcel.area;
            Fy += (sheerY2 - sheerY1) * realWorldParams.dynamicViscosity * parcel.area;
            u1x[arrayIndex(x, y)] += realToModel * Fx / mass;
            u1y[arrayIndex(x, y)] += realToModel * Fy / mass;
        }
    }
}

function draw() {
    var dataIndex, imgIndex;
    for (var y = 0; y < SIZE; y++) {
        for (var x = 0; x < SIZE; x++) {
            dataIndex = arrayIndex(x, y);
            imgIndex = dataIndex * 4;
            imageData.data[imgIndex + 0] = (visualisationMode === 1 || visualisationMode === 4) ? p0[dataIndex]   * visualScale[1] + visualOffset[1] : 0;
            imageData.data[imgIndex + 1] = (visualisationMode === 3 || visualisationMode === 4) ? u0x[dataIndex]  * visualScale[3] + visualOffset[3] : 0;
            imageData.data[imgIndex + 2] = (visualisationMode === 3 || visualisationMode === 4) ? u0y[dataIndex]  * visualScale[3] + visualOffset[3] : 0;
            imageData.data[imgIndex + 3] = (visualisationMode === 2 || visualisationMode === 4) ? div[dataIndex]  * visualScale[2] + visualOffset[2] : 255;
            // Add time lines / streaklines
        }
    }
    ctx.putImageData(imageData, 0, 0);
    replayFrames.push(ctx.getImageData(0, 0, SIZE, SIZE));
    document.getElementById("replayButton").innerHTML = "Replay " + replayFrames.length + " frames";
    reportSample();
}

function updateView(view){
    visualisationMode = view;
    for (var field = 1, maxField, minField, sampleVal, sampleField, averageVal = 0; field <= 3; field++){
        switch (field) {
            case 1: sampleField = p0    ;console.log("p0");     break;
            case 2: sampleField = div   ;console.log("div");    break;
            case 3: sampleField = u0x   ;console.log("u0y");    break;
        }
        minField =  9999999;
        maxField = -9999999;
        for (var y = 0; y < SIZE; y++) {
            for (var x = 0; x < SIZE; x++) {
                sampleVal = sampleField[arrayIndex(x, y)] || 0;
                minField = Math.min(minField, sampleVal);
                maxField = Math.max(maxField, sampleVal);
                averageVal += sampleVal;
            }
        }
        averageVal /= SIZE * SIZE;
        console.log("averageVal : " + averageVal);
        console.log("minField   : " + minField);
        console.log("maxField   : " + maxField);
        if (maxField === minField) {
            visualScale[field] = 1000;
        } else {
            visualScale[field] = 255 / (maxField - minField);
        }
        visualOffset[field] = - minField * visualScale[field];
    }
}
function contrastBoost(){
    if (visualisationMode === 4) visualScale[2] *= 2;
    else visualScale[visualisationMode] *= 2;
}
function iterateDeltaT(){
    physics();
    draw(); // This is more crude than bramble double buffer, as hard coding 0th frame
    rafRef = window.requestAnimationFrame(iterateDeltaT);
}
function replay(){
    window.cancelAnimationFrame(rafRef);
    if (replayIndex < replayFrames.length - 1){
        ctx.putImageData(replayFrames[replayIndex++], 0, 0);
        document.getElementById("replayButton").innerHTML = replayIndex + " of " + replayFrames.length;
        rafRef = window.requestAnimationFrame(replay);
    } else {
        replayIndex = 0;
        rafRef = window.requestAnimationFrame(iterateDeltaT);
    }
}
updateView(1);
