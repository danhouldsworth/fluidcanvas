/*
simpleTunnel - Navier Stokes 2D solver on CPU / Canvas.
Aims :

ADD function for returning shape boolean - Also think about how return surface, so can get pressure vector / friction vector into L / D / moment
ADD steakline view mode
ADD timeline view mode
Add vector arrow view mode

Output - Lift / Drag / Moment of object
*/

"use strict";
/* global replay */
/* jshint
browser       : true,
bitwise : false
*/

var canvas  = document.getElementById("canvas"),
    ctx     = canvas.getContext("2d"),
    mouseSample = document.getElementById("mouseSample"),
    visualisationMode   = 1,
    jacobiIterations    = 2,
    STARTSIZE           = 32,
    SIZE = canvas.width = canvas.height = STARTSIZE,
    rafRef, imageData, p0, p1, u0x, u0y, u1x, u1y, div,
    replayIndex     = 0,
    mouseX          = 0,
    mouseY          = 0,
    replayFrames    = [],
    visualScale     = [],
    visualOffset    = [],
    realWorldParams = {                         // ** All in SI
        // density                 : 1.292,        // kg/m3                        @ STP
        // kinematicViscosity      : 0.0000133,    // m2/s                         @ STP
        // dynamicViscosity        : 0.0000172,    // kg/m.s OR Pa.s OR N.s/m2     @ STP
        tunnelLength            : 5,            // m
        windSpeed               : 15,           // m/s
        soundSpeed              : 330           // m/s                          @ STP
    },
    parcel = {                              // ** All in SI **
        size    : (realWorldParams.tunnelLength / SIZE),
        area    : (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE),
        volume  : (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE)
    },
    modelParams = {                             // *** Length / Area / Volume in Parcels *** Time in iteratinons ***
        soundSpeed  : jacobiIterations,
        tunnelLength: SIZE,
        // deltaT      : 1 / (SIZE * realWorldParams.soundSpeed / realWorldParams.tunnelLength), // (Real world fluid) seconds per iteration ~ 15micro seconds
        windSpeed : jacobiIterations * realWorldParams.windSpeed / realWorldParams.soundSpeed // Pressure moves one step per jacobiIteration
    },
    trunc3dp    = function(raw) {return (raw * 1000 | 0) / 1000;},
    arrayIndex  = function(x, y){return (x + y * SIZE);},
    speedFromVector = function(vx, vy){return Math.sqrt(vx*vx + vy*vy);},
    sampleField = function(e) {
        mouseX = (e.clientX - canvas.getBoundingClientRect().left) | 0;
        mouseY = (e.clientY - canvas.getBoundingClientRect().top ) | 0;
        reportSample();
    },
    reportSample = function(){
        mouseSample.innerHTML = "<p>At mouse pointer :</p>";
        mouseSample.innerHTML += "<p>Pressure delta: " + trunc3dp(1000*p0[arrayIndex(mouseX, mouseY)]) + " <strong>mPa</strong></p>";
        mouseSample.innerHTML += "<p>Vx : " + trunc3dp(modelParams.modelToReal * u0x[arrayIndex(mouseX, mouseY)]) + " <strong>m/s</strong> (" + trunc3dp(u0x[arrayIndex(mouseX, mouseY)]) + " px/calc)</p>";
        mouseSample.innerHTML += "<p>Vy : " + trunc3dp(modelParams.modelToReal * u0y[arrayIndex(mouseX, mouseY)]) + " <strong>m/s</strong> (" + trunc3dp(u0y[arrayIndex(mouseX, mouseY)]) + " px/calc)</p>";
        mouseSample.innerHTML += "<p>Mach : " + trunc3dp(speedFromVector(u0x[arrayIndex(mouseX, mouseY)], u0y[arrayIndex(mouseX, mouseY)]) / modelParams.soundSpeed) + "</p>";
    };

canvas.addEventListener("mousemove", sampleField);
canvas.addEventListener("mousedown", sampleField);

function resizeArray(newSize){
    SIZE = newSize;
    // ** All in SI **
    parcel.size    = (realWorldParams.tunnelLength / SIZE);
    parcel.area    = (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE);
    parcel.volume  = (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE) * (realWorldParams.tunnelLength / SIZE);

    // *** Length / Area / Volume in Parcels *** Time in iteratinons ***
    modelParams.soundSpeed  = jacobiIterations;
    modelParams.tunnelLength= SIZE;
    modelParams.deltaT      = 1 / (SIZE * realWorldParams.soundSpeed / realWorldParams.tunnelLength); // (Real world fluid) seconds per iteration ~ 15micro seconds
    modelParams.windSpeed = jacobiIterations * realWorldParams.windSpeed / realWorldParams.soundSpeed; // Pressure moves one step per jacobiIteration

    realWorldParams.tunnelMach              = realWorldParams.windSpeed     / realWorldParams.soundSpeed;
    realWorldParams.soundTunnelLengthTime   = realWorldParams.tunnelLength  / realWorldParams.soundSpeed;
    realWorldParams.windTunnelLengthTime    = realWorldParams.tunnelLength  / realWorldParams.windSpeed;
    modelParams.tunnelMach                  = modelParams.windSpeed         / modelParams.soundSpeed;
    modelParams.soundTunnelLengthTime       = modelParams.tunnelLength      / modelParams.soundSpeed;
    modelParams.windTunnelLengthTime        = modelParams.tunnelLength      / modelParams.windSpeed;
    modelParams.deltaT                      = realWorldParams.soundTunnelLengthTime / modelParams.soundTunnelLengthTime;
    modelParams.realToModel                 = modelParams.deltaT            / parcel.size;         // seconds per ITER / meters per parcel ==> real speed (m/s) to parcel/iteration
    modelParams.modelToReal                 = 1 / modelParams.realToModel;

    var arraySize = SIZE * SIZE;
    var newp0  = new Float32Array(arraySize);
    var newp1  = new Float32Array(arraySize);
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    var newu0x = new Float32Array(arraySize);
    var newu0y = new Float32Array(arraySize);
    var newu1x = new Float32Array(arraySize);
    var newu1y = new Float32Array(arraySize);
    var newdiv = new Float32Array(arraySize);
    canvas.width = canvas.height = SIZE;
    ctx.fillRect(0, 0, SIZE, SIZE);
    imageData = ctx.getImageData(0, 0, SIZE, SIZE);
    // *****
    if (newSize === STARTSIZE) {
        for(var i = 0; i < arraySize; i++) {
            newp1[i]   = newp0[i]     = 0; // Define pressure as delta to ambient
            newu1x[i]  = newu0x[i]    = modelParams.windSpeed;
            newu1y[i]  = newu0y[i]    = 0;
        }
    } else {
        for(var x = 0; x < SIZE; x++) {
            for(var y = 0; y < SIZE; y++) {
                var newIndex = x + y * SIZE;
                var oldIndex = Math.floor(x/2) + Math.floor(y/2) * SIZE / 2; // Assumes we have just enhanced by *2 - this will fail going back down
                newp1[newIndex]   = newp0[newIndex]     = p0[oldIndex];
                newu1x[newIndex]  = newu0x[newIndex]    = u0x[oldIndex];
                newu1y[newIndex]  = newu0y[newIndex]    = u0y[oldIndex];
            }
        }
    }
    p0 = newp0;
    p1 = newp1;
    u0x = newu0x;
    u0y = newu0y;
    u1x = newu1x;
    u1y = newu1y;
    div = newdiv;
}


function physics(){

    tunnelBoundary(u0x, u0y);
    objectBoundary(u0x, u0y);

    advection(u0x, u1x);
    advection(u0y, u1y);
    // diffusion();
    // we now have u* or w as a intermediary divergent vector field of velocities
    computeDivergence();
    tunnelBoundary(u1x, u1y); // Mainly for divergence and pressure implications
    objectBoundary(u1x, u1y);
    jacobi(p1, p0, div, -1, 4, jacobiIterations);
    // Now we project our divergent vel field w onto divergent free by HHDT and subtract the pressure gradient
    subtractPressureGradient(p0); // An even number of iterations will leave the latest pressure field in p0
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
    // for (var y = 0; y < SIZE; y++) {
    //     for (var x = 0; x < SIZE; x++) {
    //         ux[arrayIndex(x, y)] += 0.0001; // 500 iterations to get to flow speed
    //     }
    // }

    // NoSlip / NoThru @ top and bottom sides
    for (var x = 0; x < SIZE - 1; x++) {
        // index = x + y * SIZE;
        ux[arrayIndex(x, 0)] = 0;
        uy[arrayIndex(x, 0)] = 0;
        ux[arrayIndex(x, 1)] = 0;
        uy[arrayIndex(x, 1)] = 0;
        div[arrayIndex(x, 1)] = 0;
        ux[arrayIndex(x, SIZE - 1)] = 0;
        uy[arrayIndex(x, SIZE - 1)] = 0;
        ux[arrayIndex(x, SIZE - 2)] = 0;
        uy[arrayIndex(x, SIZE - 2)] = 0;
        div[arrayIndex(x, SIZE - 2)] = 0;
    }
    // NoSlip / Laminer flow speed at left / right edges
    for (var y = 2; y < SIZE - 3; y++) {
        ux[arrayIndex(0, y)] = modelParams.windSpeed;
        uy[arrayIndex(0, y)] = 0.0;
        ux[arrayIndex(1, y)] = modelParams.windSpeed;
        uy[arrayIndex(1, y)] = 0.0;
        div[arrayIndex(1, y)] = 0.0;
        ux[arrayIndex(SIZE - 1, y)] = modelParams.windSpeed;
        uy[arrayIndex(SIZE - 1, y)] = 0.0;
        ux[arrayIndex(SIZE - 2, y)] = modelParams.windSpeed;
        uy[arrayIndex(SIZE - 2, y)] = 0.0;
        div[arrayIndex(SIZE - 2, y)] = 0.0;
    }
}

function objectBoundary(ux, uy){

    var x0  = Math.floor(SIZE / 3),
        y0  = Math.floor(SIZE / 2),
        r   = Math.floor(SIZE / 40),
        index;
    // -- CIRCLE
    for (var y = y0 - r; y < y0 + r; y++) {
        for (var x = x0 - r; x < x0 + r; x++) {
            index = x + y * SIZE;
            if ((x-x0)*(x-x0)+(y-y0)*(y-y0) < (r*r) ) {
                ux[index] = 0;
                uy[index] = 0;
                div[index] = 0;
            }
        }
    }

}

function advection(src, dest) {

    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION

    function bilerp(x, y) {
        // x,y as fractional pixel index
        // Linear interpolation between 2 values a-->b, parametised by 0<c<1
        var lerp = function (a, b, c) {return a * (1 - c) + b * c;};
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

    for (var y = 1, index; y < SIZE - 1; y++) {
        for (var x = 1, dx, dy; x < SIZE - 1; x++) {
            index = x + y * SIZE;
            dx = u0x[index];                         // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
            dy = u0y[index];                         // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
            dest[index] = bilerp(x - dx, y - dy);    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
        }
    }
}

function computeDivergence() {
    // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
    // Divergence of intermediary velocity field.
    for (var y = 1, index; y <= SIZE - 2; y++) {
        for (var x = 1; x <= SIZE - 2; x++) {
            index = x + y * SIZE;
            div[index]  = u0x[index + 1];
            div[index] -= u0x[index - 1];
            div[index] += u0y[index + SIZE];
            div[index] -= u0y[index - SIZE];
            div[index] *= 0.5; // /2*deltaX
        }
    }
}

function jacobi(dest, src, b, alpha, beta, iterations) {
    for(var i = 0; i < iterations; i++) {
        for(var y = 1, index; y <= SIZE-2; y++) {
            for(var x = 1; x <= SIZE-2; x++) {
                index = x + y * SIZE;
                dest[index]  = src[index + 1];
                dest[index] += src[index - 1];
                dest[index] += src[index + SIZE];
                dest[index] += src[index - SIZE];
                dest[index] += alpha * b[index];
                dest[index] /= beta;
            }
        }
        var temp = src;
        src = dest;
        dest = temp;
    }
}

function subtractPressureGradient(p) {
    for (var y = 1, index; y <= SIZE - 2; y++) {
        for (var x = 1; x <= SIZE - 2; x++) {
            index = x + y * SIZE;
            u1x[index] -= 0.5 * (p[index + 1]     - p[index - 1]);
            u1y[index] -= 0.5 * (p[index + SIZE]  - p[index - SIZE]);
        }
    }
}

// function diffusion(){
//     var mass        = realWorldParams.density * parcel.volume;
//     var realToModel = modelParams.deltaT / parcel.size; // seconds per ITER / meters per parcel ==> real speed (m/s) to parcel/iteration
//     // **** ALL VELOCITITES & FLUX ARE IN PARCELs / ITERATION
//     for (var y = 1, index; y < SIZE - 1; y++) {
//         for(var x = 1, sheerX0, sheerX1, sheerX2, sheerY0, sheerY1, sheerY2, Fx, Fy; x < SIZE - 1; x++) {
//             index = x + y * SIZE;
//             sheerY0 = u0y[arrayIndex(x - 1, y)];
//             sheerY1 = u0y[arrayIndex(x,     y)];
//             sheerY2 = u0y[arrayIndex(x + 1, y)];
//             sheerX0 = u0x[arrayIndex(x,     y - 1)];
//             sheerX1 = u0x[arrayIndex(x,     y)];
//             sheerX2 = u0x[arrayIndex(x,     y + 1)];
//             // dynamicViscosity applies Force on fluid parcel in direction of strain force
//             Fx  = (sheerX0 - sheerX1) * realWorldParams.dynamicViscosity * parcel.area;
//             Fx += (sheerX2 - sheerX1) * realWorldParams.dynamicViscosity * parcel.area;
//             Fy  = (sheerY0 - sheerY1) * realWorldParams.dynamicViscosity * parcel.area;
//             Fy += (sheerY2 - sheerY1) * realWorldParams.dynamicViscosity * parcel.area;
//             u1x[index] += realToModel * Fx / mass;
//             u1y[index] += realToModel * Fy / mass;
//         }
//     }
// }

function draw() {
    var dataIndex, imgIndex;
    for (var y = 0; y <= SIZE-1; y++) {
        for (var x = 0; x <= SIZE-1; x++) {
            dataIndex = x + y * SIZE;
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
        for (var y = 0; y <= SIZE-1; y++) {
            for (var x = 0; x <= SIZE-1; x++) {
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
    if (visualisationMode === 4) {visualScale[2] *= 2;} // Dial up flux if rendering composite of all
    else {visualScale[visualisationMode] *= 2;}
}
function iterateDeltaT(){
    physics();physics();physics();physics();physics();physics();physics();physics();physics();physics(); // 10 CalcsPerRender
    draw(); // This is more crude than bramble double buffer, as hard coding 0th frame
    rafRef = window.requestAnimationFrame(iterateDeltaT);
}
function iterationsChange(delta){
    jacobiIterations += delta;
    document.getElementById("iterations").innerHTML = "Jacobie iterations = " + jacobiIterations;
}
function resolutionChange(delta){
    resizeArray(SIZE * delta);
    document.getElementById("resolution").innerHTML = "Enhance this! (" + SIZE + "x" + SIZE + ")";
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

resolutionChange(1);
updateView(1);iterationsChange(0);