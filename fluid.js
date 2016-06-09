var canvas  = document.getElementById('canvas'),
    ctx     = canvas.getContext('2d'),
    WIDTH   = canvas.width  = 200,
    HEIGHT  = canvas.height = 200,
    sx      = canvas.width / canvas.clientWidth,
    sy      = canvas.height / canvas.clientHeight,
    offsetLeft  = canvas.getBoundingClientRect().left,
    offsetTop   = canvas.getBoundingClientRect().top,
    mouseX      = 0,
    mouseY      = 0,
    lastMouseX  = mouseX,
    lastMouseY  = mouseY,
    imageData;

canvas.addEventListener('mousemove', function(e) {
    mouseX = (e.clientX - offsetLeft)   | 0;
    mouseY = (e.clientY - offsetTop)    | 0;
});

ctx.fillRect(0, 0, WIDTH, HEIGHT);
imageData = ctx.getImageData(0, 0, WIDTH, HEIGHT);

var velocityField0  = new Float32Array(WIDTH * HEIGHT * 2),
    velocityField1  = new Float32Array(WIDTH * HEIGHT * 2),
    pressureField0  = new Float32Array(WIDTH * HEIGHT),
    pressureField1  = new Float32Array(WIDTH * HEIGHT),
    divergenceField = new Float32Array(WIDTH * HEIGHT),
    u0x = sampler(velocityField0,   2, 0),
    u0y = sampler(velocityField0,   2, 1),
    u1x = sampler(velocityField1,   2, 0),
    u1y = sampler(velocityField1,   2, 1),
    p0  = sampler(pressureField0,   1, 0),
    p1  = sampler(pressureField1,   1, 0),
    div = sampler(divergenceField,  1, 0),
    step = 4.0;

for (var i = 0; i < pressureField0.length; i++) {
    pressureField1[i] = pressureField0[i] = 0;
}
for(i = 0; i < velocityField0.length; i++) {
    velocityField1[i] = velocityField0[i] = (Math.random()-0.5)*1.0;
}
// velocityboundary(u0x, u0y);

function simulate(){

    // velocityboundary(u0x, u0y);
    advect(u0x, u0y, u0x, u1x, step);
    advect(u0x, u0y, u0y, u1y, step);
    addMouseForce(u1x, u1y);
    computeDivergence(u1x, u1y, div);
    // needs an even number of iterations
    fastjacobi(p0, p1, div, 16);
    subtractPressureGradient(u1x, u1y, p0);

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

function addMouseForce(ux, uy) {
    var x   = mouseX * sx,
        y   = mouseY * sy,
        dx  = mouseX - lastMouseX,
        dy  = mouseY - lastMouseY;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
    ux(x, y, ux(x, y) - dx * 2);
    uy(x, y, uy(x, y) - dy * 2);
}

function velocityboundary(ux, uy) {
    for (var x = 0; x < WIDTH - 1; x++) {
        ux(x, 0, -ux(x, 1));
        uy(x, 0, -uy(x, 1));

        ux(x, HEIGHT - 1, -ux(x, HEIGHT - 2));
        uy(x, HEIGHT - 1, -uy(x, HEIGHT - 2));
    }
    for (var y = 0; y < HEIGHT - 1; y++) {
        ux(0, y, -ux(1, y));
        uy(0, y, -uy(1, y));

        ux(WIDTH - 1, y, -ux(WIDTH - 2, y));
        uy(WIDTH - 1, y, -uy(WIDTH - 2, y));
    }
}

function lerp(a, b, c) {
    c = c < 0 ? 0 : (c > 1 ? 1 : c);
    return a * (1 - c) + b * c;
}

function sampler(field, stride, offset) {
    var f = function(x, y, value) {
        x = (x < 0 ? 0 : (x > WIDTH  - 1 ? WIDTH  - 1 : x)) | 0;
        y = (y < 0 ? 0 : (y > HEIGHT - 1 ? HEIGHT - 1 : y)) | 0;
        var index = (x + y * WIDTH) * stride + offset;
        if (value)  field[index] = value;
        else        return field[index];
    };
    f.field = field;
    return f;
}

function bilerp(sample, x, y) {
    var x0 = ~~x,
        y0 = ~~y,
        x1 = x0 + 1,
        y1 = y0 + 1,
        p00 = sample(x0, y0),
        p01 = sample(x0, y1),
        p10 = sample(x1, y0),
        p11 = sample(x1, y1);
    return lerp(lerp(p00, p10, x - x0), lerp(p01, p11, x - x0), y - y0);
}

function advect(ux, uy, src, dest, t) {
    for (var y = 1; y < HEIGHT - 2; y++) {
        for (var x = 1; x < WIDTH - 2; x++) {
            var vx = ux(x, y) * t,
                vy = uy(x, y) * t;
            dest(x, y, bilerp(src, x + vx, y + vy));
        }
    }
}

function computeDivergence(ux, uy, div) {
    for (var y = 1; y < HEIGHT - 2; y++) {
        for (var x = 1, index; x < WIDTH - 2; x++) {
            // compute divergence using central difference
            index = x + y * WIDTH;
            div.field[index]  = ux(x + 1,  y);
            div.field[index] -= ux(x - 1,  y);
            div.field[index] += uy(x,      y + 1);
            div.field[index] -= uy(x,      y - 1);
            div.field[index] /= 2;
        }
    }
}

function fastjacobi(p0, p1, div, iterations) {
    for (i = 0; i < iterations; i++) {
        for (var y = 1; y < HEIGHT - 2; y++) {
            for (var x = 1, index; x < WIDTH - 2; x++) {
                index = x + y * WIDTH;
                p1.field[index] = -div.field[index];
                p1.field[index] += p0.field[index - 1];
                p1.field[index] += p0.field[index + 1];
                p1.field[index] += p0.field[index - WIDTH];
                p1.field[index] += p0.field[index + WIDTH];
                p1.field[index] /= 4;
            }
        }
        var temp = p0;
        p0 = p1;
        p1 = temp;
    }
}

function subtractPressureGradient(ux, uy, p) {
    for (var y = 1; y < HEIGHT - 2; y++) {
        for(var x = 1, x0, x1, y0, y1, dx, dy; x < WIDTH - 2; x++) {
            x0 = p(x - 1,   y);
            x1 = p(x + 1,   y);
            y0 = p(x,       y - 1);
            y1 = p(x,       y + 1);
            dx = (x1 - x0) / 2;
            dy = (y1 - y0) / 2;

            ux(x, y, ux(x, y) - dx);
            uy(x, y, uy(x, y) - dy);
        }
    }
}

function draw(ux, uy, p) {
    for (var y = 0; y < HEIGHT - 1; y++) {
        for (var x = 0, dataIndex, pIndex, uIndex; x < WIDTH - 1; x++) {
            pIndex = (y * WIDTH + x);
            uIndex = pIndex * 2;
            dataIndex = pIndex * 4;
            // imageData.data[dataIndex + 0] = div(x, y) * 128 + 128;
            imageData.data[dataIndex + 1] = ux(x, y) * 128 + 128;
            imageData.data[dataIndex + 2] = uy(x, y) * 128 + 128;
            // imageData.data[dataIndex + 3] = p(x, y)  * 255;
        }
    }
    ctx.putImageData(imageData, 0, 0);
}

function animate(){
    simulate();
    draw(u0x, u0y, p0);
    window.requestAnimationFrame(animate);
};

animate();