function generateEllipsoidMesh(centerX, centerY, centerZ, radX, radY, radZ, subdivisions) {
  var vertices = [];
  var indices = [];
  var normals = [];

  // Loop untuk sudut theta (elevasi)
  for (var i = 0; i <= subdivisions; i++) {
      var theta = i * Math.PI / subdivisions;  // Dari 0 hingga PI
      var sinTheta = Math.sin(theta);
      var cosTheta = Math.cos(theta);

      // Loop untuk sudut phi (azimut)
      for (var j = 0; j <= subdivisions; j++) {
          var phi = j * 2 * Math.PI / subdivisions;  // Dari 0 hingga 2*PI
          var sinPhi = Math.sin(phi);
          var cosPhi = Math.cos(phi);

          // Hitung koordinat x, y, dan z untuk titik pada permukaan ellipsoid
          var x = centerX + radX * sinTheta * cosPhi;  // Sumbu x
          var y = centerY + radY * cosTheta;  // Sumbu y
          var z = centerZ + radZ * sinTheta * sinPhi;  // Sumbu z

          // Hitung normal titik
          var nx = (x - centerX) / radX;
          var ny = (y - centerY) / radY;
          var nz = (z - centerZ) / radZ;

          // Normalisasi normal
          var normLength = Math.sqrt(nx * nx + ny * ny + nz * nz);
          nx /= normLength;
          ny /= normLength;
          nz /= normLength;

          // Tambahkan titik dan normal ke array
          vertices.push(x, y, z);
          normals.push(nx, ny, nz);
      }
  }

  // Buat indeks untuk menggambar segitiga
  for (var i = 0; i < subdivisions; i++) {
      for (var j = 0; j < subdivisions; j++) {
          var first = i * (subdivisions + 1) + j;
          var second = first + subdivisions + 1;

          // Segitiga pertama
          indices.push(first);
          indices.push(second);
          indices.push(first + 1);

          // Segitiga kedua
          indices.push(second);
          indices.push(second + 1);
          indices.push(first + 1);
      }
  }

  return {
      vertices: vertices,
      indices: indices,
      normals: normals
  };
}

function generateCone(x, y, z, rad, height, segments, color = [1, 1, 1], isRainbow = false) {
  var vertices = [];
  var colors = [];

  // Center point
  vertices.push(x, y, z);

  colors.push(color[0], color[1], color[2]);

  if(isRainbow){
    var rainbow = [
      [1.0,0.0,1.0],
      [0.0,1.0,0.0],
      [1.0,0.0,1.0]
    ]
  }

  // Vertices around the base
  for (var i = 0; i < segments; i++) {
      var angle = (i / segments) * Math.PI * 2;
      var vx = rad * Math.cos(angle) + x;
      var vy = rad * Math.sin(angle) + y;
      var vz = 0;
      vertices.push(vx, vy, vz);
      colors.push(color[0], color[1], color[2]);

      if(isRainbow){
        var colorIndex = i % rainbow.length;
        colors = colors.concat(rainbow[colorIndex]);
      }

  }
  
  // Vertex at the apex
  vertices.push(x, y, height);
  colors.push(color[0], color[1], color[2]);

  // Faces
  var faces = [];
  for (var i = 1; i <= segments; i++) {
      faces.push(0, i, (i % segments) + 1);
      faces.push(segments + 1, (i % segments) + 1, i);
      faces.push(i, (i % segments) + 1, segments + 1);
  }
  
  return { vertices: vertices, colors: colors, faces: faces };
}
function generateCone2(x, z, y, rad, height, segments, color = [1, 1, 1], isRainbow = false) {
  var vertices = [];
  var colors = [];

  // Center point
  vertices.push(x, y, z);

  colors.push(color[0], color[1], color[2]);

  if(isRainbow){
    var rainbow = [
      [1.0,0.0,1.0],
      [0.0,1.0,0.0],
      [1.0,0.0,1.0]
    ]
  }

  // Vertices around the base
  for (var i = 0; i < segments; i++) {
      var angle = (i / segments) * Math.PI * 2;
      var vx = rad * Math.cos(angle) + x;
      var vy = rad * Math.sin(angle) + y;
      var vz = 0;
      vertices.push(vx, vy, vz);
      colors.push(color[0], color[1], color[2]);

      if(isRainbow){
        var colorIndex = i % rainbow.length;
        colors = colors.concat(rainbow[colorIndex]);
      }

  }
  
  // Vertex at the apex
  vertices.push(x, y, height);
  colors.push(color[0], color[1], color[2]);

  // Faces
  var faces = [];
  for (var i = 1; i <= segments; i++) {
      faces.push(0, i, (i % segments) + 1);
      faces.push(segments + 1, (i % segments) + 1, i);
      faces.push(i, (i % segments) + 1, segments + 1);
  }
  
  return { vertices: vertices, colors: colors, faces: faces };
}

var Icosahedron3D = (function () {
  function Icosahedron3D(quality) {
      this._quality = quality;
      this._calculateGeometry();
  }
  Icosahedron3D.prototype._calculateGeometry = function () {
      this.Points = [];
      this.TriangleIndices = [];
      this._middlePointIndexCache = {};
      this._index = 0;
      var t = (1.0 + Math.sqrt(5.0)) / 2.0;
      this._addVertex(-1, t, 0);
      this._addVertex(1, t, 0);
      this._addVertex(-1, -t, 0);
      this._addVertex(1, -t, 0);
      this._addVertex(0, -1, t);
      this._addVertex(0, 1, t);
      this._addVertex(0, -1, -t);
      this._addVertex(0, 1, -t);
      this._addVertex(t, 0, -1);
      this._addVertex(t, 0, 1);
      this._addVertex(-t, 0, -1);
      this._addVertex(-t, 0, 1);
      this._addFace(0, 11, 5);
      this._addFace(0, 5, 1);
      this._addFace(0, 1, 7);
      this._addFace(0, 7, 10);
      this._addFace(0, 10, 11);
      this._addFace(1, 5, 9);
      this._addFace(5, 11, 4);
      this._addFace(11, 10, 2);
      this._addFace(10, 7, 6);
      this._addFace(7, 1, 8);
      this._addFace(3, 9, 4);
      this._addFace(3, 4, 2);
      this._addFace(3, 2, 6);
      this._addFace(3, 6, 8);
      this._addFace(3, 8, 9);
      this._addFace(4, 9, 5);
      this._addFace(2, 4, 11);
      this._addFace(6, 2, 10);
      this._addFace(8, 6, 7);
      this._addFace(9, 8, 1);
      this._refineVertices();
  };
  Icosahedron3D.prototype._addVertex = function (x, y, z) {
      var length = Math.sqrt(x * x + y * y + z * z);
      this.Points.push({
          x: x / length,
          y: y / length,
          z: z / length
      });
      return this._index++;
  };
  Icosahedron3D.prototype._addFace = function (x, y, z) {
      this.TriangleIndices.push(x);
      this.TriangleIndices.push(y);
      this.TriangleIndices.push(z);
  };
  Icosahedron3D.prototype._refineVertices = function () {
      for (var i = 0; i < this._quality; i++) {
          var faceCount = this.TriangleIndices.length;
          for (var face = 0; face < faceCount; face += 3) {
              var x1 = this.TriangleIndices[face];
              var y1 = this.TriangleIndices[face + 1];
              var z1 = this.TriangleIndices[face + 2];
              var x2 = this._getMiddlePoint(x1, y1);
              var y2 = this._getMiddlePoint(y1, z1);
              var z2 = this._getMiddlePoint(z1, x1);
              this._addFace(x1, x2, z2);
              this._addFace(y1, y2, x2);
              this._addFace(z1, z2, y2);
              this._addFace(x2, y2, z2);
          }
      }
  };
  Icosahedron3D.prototype._getMiddlePoint = function (p1, p2) {
      var firstIsSmaller = p1 < p2;
      var smallerIndex = firstIsSmaller ? p1 : p2;
      var greaterIndex = firstIsSmaller ? p2 : p1;
      var key = (smallerIndex << 32) + greaterIndex;
      var p = this._middlePointIndexCache[key];
      if (p !== undefined)
          p;
      var point1 = this.Points[p1];
      var point2 = this.Points[p2];
      var middle = {
          x: (point1.x + point2.x) / 2.0,
          y: (point1.y + point2.y) / 2.0,
          z: (point1.z + point2.z) / 2.0,
      };
      var i = this._addVertex(middle.x, middle.y, middle.z);
      this._middlePointIndexCache[key] = i;
      return i;
  };
  return Icosahedron3D;
})();

var Cylinder3D = (function () {
  function Cylinder3D(radius, height, segmentsHeight, segmentsRadius, closed = true) {
    this._radius = radius;
    this._height = height;
    this._segmentsHeight = segmentsHeight;
    this._segmentsRadius = segmentsRadius;
    this._closed = closed;
    this._calculateGeometry();
  }

  Cylinder3D.prototype._calculateGeometry = function () {
    this.Points = [];
    this.TriangleIndices = [];

    var halfHeight = this._height / 2.0;
    var angleStep = (Math.PI * 2) / this._segmentsRadius;
    var segmentHeightStep = this._height / this._segmentsHeight;

    // Generate vertices on the top and bottom circles
    for (var i = 0; i <= this._segmentsRadius; i++) {
      var angle = i * angleStep;
      var x = this._radius * Math.cos(angle);
      var z = this._radius * Math.sin(angle);

      this.Points.push({ x: x, y: halfHeight, z: z }); // Top circle
      this.Points.push({ x: x, y: -halfHeight, z: z }); // Bottom circle
    }

    // Generate triangles for the body
    for (var i = 0; i < this._segmentsRadius; i++) {
      for (var j = 0; j < this._segmentsHeight; j++) {
        var topIndex1 = i * 2;
        var topIndex2 = (i + 1) * 2;
        var bottomIndex1 = topIndex1 + 1;
        var bottomIndex2 = topIndex2 + 1;

        // Connect top and bottom quadrilaterals of each segment
        this.TriangleIndices.push(topIndex1, bottomIndex1, topIndex2);
        this.TriangleIndices.push(topIndex2, bottomIndex1, bottomIndex2);
      }
    }

    // Add top and bottom caps (if closed)
    if (this._closed) {
      var centerIndex = this.Points.length; // Index of the center point
      var topCenter = { x: 0, y: halfHeight, z: 0 };
      var bottomCenter = { x: 0, y: -halfHeight, z: 0 };

      // Add top and bottom center points
      this.Points.push(topCenter);
      this.Points.push(bottomCenter);

      // Generate top cap triangles
      for (var i = 0; i < this._segmentsRadius; i++) {
        var topIndex = i * 2;
        this.TriangleIndices.push(centerIndex, topIndex, (i + 1) * 2);
      }

      // Generate bottom cap triangles
      for (var i = 0; i < this._segmentsRadius; i++) {
        var bottomIndex = this.Points.length - 2 - i;
        this.TriangleIndices.push(centerIndex, bottomIndex, bottomIndex - 1);
      }
    }
  };

  return Cylinder3D;
})();


function generateHead(x, y, z, radius, eyeRadius, eyeXOffset, eyeYOffset, eyeZOffset) {
  var head = generateSphere(x, y, z, radius, 16, 16); // Generate head sphere

  // Generate eyes (spheres positioned at offsets)
  var leftEye = generateSphere(x - eyeXOffset, y + eyeYOffset, z + eyeZOffset, eyeRadius, 1000, 1000);
  var rightEye = generateSphere(x + eyeXOffset, y + eyeYOffset, z + eyeZOffset, eyeRadius, 1000, 1000);

  // Combine head, eyes, and additional features (modify as needed)
  var headMesh = {
    vertices: head.vertices.concat(leftEye.vertices, rightEye.vertices),
    colors: head.colors.concat(leftEye.colors, rightEye.colors),
    faces: head.faces.concat(leftEye.faces, rightEye.faces)
  };

  // Add a simple nose cylinder (modify for desired shape)
  var noseRadius = radius * 0.1;
  var noseHeight = radius * 0.2;
  var noseCylinder = generateCylinder(x, y + radius * 0.3, noseRadius, noseHeight, 16, [1, 0.7, 0.5]); // Adjust color as needed
  headMesh.vertices = headMesh.vertices.concat(noseCylinder.vertices);
  headMesh.colors = headMesh.colors.concat(noseCylinder.colors);
  headMesh.faces = headMesh.faces.concat(noseCylinder.faces);

  return headMesh;
}

// Function for generating a sphere (omitted for brevity)
function generateSphere(x, y, z, radius, segmentsH, segmentsV) {
  // Arrays to store data
  var vertices = [];
  var colors = [];
  var faces = [];

  // Calculate latitude (phi) and longitude (theta) angles
  for (var phi = 0; phi <= segmentsV; phi++) {
    var latitudeAngle = Math.PI * phi / segmentsV;
    var yValue = Math.sin(latitudeAngle);

    for (var theta = 0; theta <= segmentsH; theta++) {
      var longitudeAngle = Math.PI * 2 * theta / segmentsH;
      var xValue = Math.cos(longitudeAngle) * Math.sin(latitudeAngle);
      var zValue = Math.cos(latitudeAngle);

      // Calculate vertex position
      var vx = xValue * radius + x;
      var vy = yValue * radius + y;
      var vz = zValue * radius + z;

      // Add vertex to array
      vertices.push(vx, vy, vz);

      // Calculate vertex color (simple RGB based on latitude and longitude)
      var red = 0.5 + 0.5 * xValue;
      var green = 0.5 + 0.5 * yValue;
      var blue = 0.5 + 0.5 * zValue;
      colors.push(red, green, blue);

      // Add faces for sphere (triangle strips)
      if (phi < segmentsV && theta < segmentsH) {
        var face1 = [phi * segmentsH + theta, (phi + 1) * segmentsH + theta, (phi + 1) * segmentsH + (theta + 1)];
        var face2 = [phi * segmentsH + theta, (phi + 1) * segmentsH + (theta + 1), phi * segmentsH + (theta + 1)];
        faces.push(face1, face2);
      }
    }
  }

  // Return object containing generated data
  return { vertices: vertices, colors: colors, faces: faces };
}



function generateCylinder(x, y, rad, height, segments, color = [1, 1, 1], isRainbow = false) {
  var vertices = [];
  var colors = [];
  var faces = [];

  if (isRainbow) {
    var rainbow = [
      [1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0],
      [0.0, 0.0, 1.0]
    ];
  }

  // Vertices around the base (adjusted for center position)
  for (var i = 0; i < segments; i++) {
    var angle = (i / segments) * Math.PI * 2;
    var vx = rad * Math.cos(angle) + x;
    var vy = rad * Math.sin(angle) + y;
    vertices.push(vx, vy, 0); // Base at z = 0
    colors.push(color[0], color[1], color[2]);

    if (isRainbow) {
      var colorIndex = i % rainbow.length;
      colors = colors.concat(rainbow[colorIndex]);
    }
  }

  // Vertices around the top (adjusted for center position and tilt)
  for (var i = 0; i < segments; i++) {
    var angle = (i / segments) * Math.PI * 2;
    var vx = rad * Math.cos(angle) + x;
    var vy = rad * Math.sin(angle) + y; // Typo fix: sin instead of sin

    // Approximate tilt correction (adjust factor based on desired accuracy)
    var tiltFactor = 0.001; // Experiment with this value
    var vz = height + (tiltFactor * Math.abs(vx - x) + tiltFactor * Math.abs(vy - y));

    vertices.push(vx, vy, vz);
    colors.push(color[0], color[1], color[2]);

    if (isRainbow) {
      var colorIndex = i % rainbow.length;
      colors = colors.concat(rainbow[colorIndex]);
    }
  }

  // Closing face for the base (using triangle fan for better alignment)
  for (var i = 0; i < segments; i++) {
    faces.push(segments, (i + 1) % segments, i); // Center vertex, then counter-clockwise order
    colors.push(color[0], color[1], color[2]);

    if (isRainbow) {
      var colorIndex = i % rainbow.length;
      colors = colors.concat(rainbow[colorIndex]);
    }
  }

  // Faces for the side (using triangle strips for better performance)
  for (var i = 0; i < segments; i++) {
    var nextIndex = (i + 1) % segments;
    faces.push(i, nextIndex, segments + i);
    faces.push(segments + i, segments + nextIndex, nextIndex);
  }

  return { vertices: vertices, colors: colors, faces: faces };
}

function generateCube(x, y, z, size, colors) {
  var vertices = [];
  var colorArray = [];

  // Generate vertices and colors for each face of the cube
  for (var dx = -1; dx <= 1; dx += 2) {
    for (var dy = -1; dy <= 1; dy += 2) {
      for (var dz = -1; dz <= 1; dz += 2) {
        vertices.push(x + dx * size[0] / 2, y + dy * size[1] / 2, z + dz * size[2] / 2);
        colorArray.push(colors[0], colors[1], colors[2]);
      }
    }
  }

  // Define the faces of the cube (each face is a square made up of two triangles)
  var faces = [
    0, 2, 1, 1, 2, 3, // front
    4, 5, 6, 6, 5, 7, // back
    0, 1, 4, 4, 1, 5, // left
    2, 6, 3, 3, 6, 7, // right
    0, 4, 2, 2, 4, 6, // top
    1, 3, 5, 5, 3, 7  // bottom
  ];

  return { vertices: vertices, colors: colorArray, faces: faces };
}


  function generateSabitKontol(x, y, rad, height, segments, color = [1,1,1]) {
    var vertices = [];
    var colors = [];
  
    // Center point
    vertices.push(x, y, 0);
   console.log(color)
    colors.push(0,0,0);
  
    var rainbow = [
      [1.0,1.0,1.0],
      [1.0,1.0,1.0],
      [1.0,1.0,1.0]
    ]
  
    // Vertices around the base
    for (var i = 0; i < segments; i++) {
        var angle = (i / segments) * Math.PI * 2;
        var vx = rad * Math.cos(angle) + x;
        var vy = rad * Math.sin(angle) + y;
        var vz = 0;
        vertices.push(vx, vy, vz);
        colors.push(0,0,0);
  
  
        // var colorIndex = i % rainbow.length;
        // colors = colors.concat(rainbow[colorIndex]);
  
    }
  
    // Vertex at the apex
    vertices.push(x, y, height);
    colors.push(0,0,0);
  
    // Faces
    var faces = [];
    for (var i = 1; i <= segments; i++) {
        faces.push(0, i, (i % segments) + 1);
        faces.push(segments + 1, (i % segments) + 1, i);
        faces.push(i, (i % segments) + 1, segments + 1);
    }
  
    return { vertices: vertices, colors: colors, faces: faces };
  }

  function generateQuadraticBezierCurve(controlPoints, t) {
    const vertices = [];
    const faces = [];
  
    // Expected number of vertices (rounded up)
    const expectedVertices = Math.ceil(controlPoints.length * 2 / t);
  
    // Generate vertices for each point on the curve
    for (var i = 0; i <= 1; i += t) {
      const point = calculateQuadraticBezierPoint(i, controlPoints);
      vertices.push(point[0], point[1]);
    }
  
    // Truncate vertices array if it has more elements than expected
    vertices.length = Math.min(vertices.length, expectedVertices * 2);
  
    // Define faces for the curve (each face is a triangle made up of three vertices)
    for (var i = 0; i < vertices.length / 2; i++) {
      faces.push(i, i + 1, i + (i % 2 === 0 ? 2 : 1));  // Check for even index
    }
  
    return { vertices: vertices, faces: faces };
  }
  
  function calculateQuadraticBezierPoint(t, controlPoints) {
    const [p0, p1, p2] = controlPoints;
  
    const x = (1 - t) * (1 - t) * p0[0] + 2 * t * (1 - t) * p1[0] + t * t * p2[0];
    const y = (1 - t) * (1 - t) * p0[1] + 2 * t * (1 - t) * p1[1] + t * t * p2[1];
  
    return [x, y];
  }
  
  


  