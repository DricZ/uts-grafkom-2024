class Vec3 extends Float32Array{
    constructor(ini){
        super(3);
        if(ini instanceof Vec3){
            this[0] = ini[0]; this[1] = ini[1]; this[2] = ini[2];
        }else{
            this[0] = this[1] = this[2] = ini || 0;
        }
        this.isModified = true;
    }

    //----------------------------------------------
    //region XYZ Setters
        set(x,y,z){ this[0] = x; this[1] = y; this[2] = z; this.isModified = true; return this;}

        get x(){ return this[0]; }	set x(val){ this[0] = val; this.isModified = true; }
        get y(){ return this[1]; }	set y(val){ this[1] = val; this.isModified = true; }
        get z(){ return this[2]; }	set z(val){ this[2] = val; this.isModified = true; }
    //endregion

    //----------------------------------------------
    //region Methods
        magnitude(v){
            //Only get the magnitude of this vector
            if(v === undefined) return Math.sqrt( this[0]*this[0] + this[1]*this[1] + this[2]*this[2] );

            //Get magnitude based on another vector
            var x = v[0] - this[0],
                y = v[1] - this[1],
                z = v[2] - this[2];

            return Math.sqrt( x*x + y*y + z*z );
        }

        normalize(){
            var mag = Math.sqrt( this[0]*this[0] + this[1]*this[1] + this[2]*this[2] );
            this[0] /= mag;
            this[1] /= mag;
            this[2] /= mag;
            this.isModified = true;
            return this;
        }

        multi(v){
            this[0] *= v;
            this[1] *= v;
            this[2] *= v;
            this.isModified = true;
            return this;
        }

        add(v){
            this[0] += v[0];
            this[1] += v[1];
            this[2] += v[2];
            this.isModified = true;
            return this;
        }

        clone(){ return new Vec3().set(this.x,this.y,this.z); }
        
        copy(v){
            this[0] = v[0]; this[1] = v[1]; this[2] = v[2];
            this.isModified = true;
            return this;
        }
    //endregion
}

class Quaternion extends Float32Array{
    constructor(){
        super(4);
        this[0] = this[1] = this[2] = 0;
        this[3] = 1;
        this.isModified = false;
    }
    //http://in2gpu.com/2016/03/14/opengl-fps-camera-quaternion/
    //----------------------------------------------
    //region Setter/Getters
        reset(){ this[0] = this[1] = this[2] = 0; this[3] = 1; this.isModified = false; return this; }

        rx(rad){ Quaternion.rotateX(this,this,rad); this.isModified = true; return this; }
        ry(rad){ Quaternion.rotateY(this,this,rad); this.isModified = true; return this; }
        rz(rad){ Quaternion.rotateZ(this,this,rad); this.isModified = true; return this; }
        
        //ex(deg){ Quaternion.rotateX(this,this,deg * DEG2RAD); this.isModified = true; return this; }
        //ey(deg){ Quaternion.rotateY(this,this,deg * DEG2RAD); this.isModified = true; return this; }
        //ez(deg){ Quaternion.rotateZ(this,this,deg * DEG2RAD); this.isModified = true; return this; }
    //endregion

    //----------------------------------------------
    //region Static Methods
        static multi(out,a,b){
            var ax = a[0], ay = a[1], az = a[2], aw = a[3],
            bx = b[0], by = b[1], bz = b[2], bw = b[3];

            out[0] = ax * bw + aw * bx + ay * bz - az * by;
            out[1] = ay * bw + aw * by + az * bx - ax * bz;
            out[2] = az * bw + aw * bz + ax * by - ay * bx;
            out[3] = aw * bw - ax * bx - ay * by - az * bz;
            return out;
        }

        static multiVec3(out,q,v){
            var ax = q[0], ay = q[1], az = q[2], aw = q[3],
                bx = v[0], by = v[1], bz = v[2];

            out[0] = ax + aw * bx + ay * bz - az * by;
            out[1] = ay + aw * by + az * bx - ax * bz;
            out[2] = az + aw * bz + ax * by - ay * bx;
            return out;
        }

        //https://github.com/toji/gl-matrix/blob/master/src/gl-matrix/quat.js
        static rotateX(out, a, rad){
            rad *= 0.5; 

            var ax = a[0], ay = a[1], az = a[2], aw = a[3],
                bx = Math.sin(rad), bw = Math.cos(rad);

            out[0] = ax * bw + aw * bx;
            out[1] = ay * bw + az * bx;
            out[2] = az * bw - ay * bx;
            out[3] = aw * bw - ax * bx;
            return out;
        }

        static rotateY(out, a, rad) {
            rad *= 0.5; 

            var ax = a[0], ay = a[1], az = a[2], aw = a[3],
            by = Math.sin(rad), bw = Math.cos(rad);

            out[0] = ax * bw - az * by;
            out[1] = ay * bw + aw * by;
            out[2] = az * bw + ax * by;
            out[3] = aw * bw - ay * by;
            return out;
        }

        static rotateZ(out, a, rad){
            rad *= 0.5; 

            var ax = a[0], ay = a[1], az = a[2], aw = a[3],
            bz = Math.sin(rad), bw = Math.cos(rad);

            out[0] = ax * bw + ay * bz;
            out[1] = ay * bw - ax * bz;
            out[2] = az * bw + aw * bz;
            out[3] = aw * bw - az * bz;
            return out;
        }

        //https://github.com/mrdoob/three.js/blob/dev/src/math/Quaternion.js
        static setFromEuler(out,x,y,z,order){
            var c1 = Math.cos(x/2),
                c2 = Math.cos(y/2),
                c3 = Math.cos(z/2),
                s1 = Math.sin(x/2),
                s2 = Math.sin(y/2),
                s3 = Math.sin(z/2);

            switch(order){
                case 'XYZ':			
                    out[0] = s1 * c2 * c3 + c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 - s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 + s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 - s1 * s2 * s3;
                    break;
                case 'YXZ':
                    out[0] = s1 * c2 * c3 + c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 - s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 - s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 + s1 * s2 * s3;
                    break;
                case 'ZXY':
                    out[0] = s1 * c2 * c3 - c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 + s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 + s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 - s1 * s2 * s3;
                    break;
                case 'ZYX':
                    out[0] = s1 * c2 * c3 - c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 + s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 - s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 + s1 * s2 * s3;
                    break;
                case 'YZX':
                    out[0] = s1 * c2 * c3 + c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 + s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 - s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 - s1 * s2 * s3;
                    break;
                case 'XZY':
                    out[0] = s1 * c2 * c3 - c1 * s2 * s3;
                    out[1] = c1 * s2 * c3 - s1 * c2 * s3;
                    out[2] = c1 * c2 * s3 + s1 * s2 * c3;
                    out[3] = c1 * c2 * c3 + s1 * s2 * s3;
                    break;
            }
        }
    //endregion
}

class Matrix4 extends Float32Array{
    constructor(){ super(16); this[0] = this[5] = this[10] = this[15] = 1; }  //Setup Identity

    //----------------------------------------------
    //region Methods
        translate(ary){	Matrix4.translate(this,ary[0],ary[1],ary[2]); return this;}
        resetTranslation(){ this[12] = this[13] = this[14] = 0; this[15] = 1; return this; }

        //reset data back to identity.
        reset(){ 
            for(var i=0; i <= this.length; i++) this[i] = (i % 5 == 0)? 1 : 0; //only positions 0,5,10,15 need to be 1 else 0
            return this;
        }
    //endregion

    //----------------------------------------------
    //region Static
        static identity(out){
            for(var i=0; i <= out.length; i++) out[i] = (i % 5 == 0)? 1 : 0; //only positions 0,5,10,15 need to be 1 else 0
        }

        static perspective(out, fovy, aspect, near, far){
            var f = 1.0 / Math.tan(fovy / 2),
                nf = 1 / (near - far);
            out[0] = f / aspect;
            out[1] = 0;
            out[2] = 0;
            out[3] = 0;
            out[4] = 0;
            out[5] = f;
            out[6] = 0;
            out[7] = 0;
            out[8] = 0;
            out[9] = 0;
            out[10] = (far + near) * nf;
            out[11] = -1;
            out[12] = 0;
            out[13] = 0;
            out[14] = (2 * far * near) * nf;
            out[15] = 0;
        }

        static ortho(out, left, right, bottom, top, near, far) {
            var lr = 1 / (left - right),
                bt = 1 / (bottom - top),
                nf = 1 / (near - far);
            out[0] = -2 * lr;
            out[1] = 0;
            out[2] = 0;
            out[3] = 0;
            out[4] = 0;
            out[5] = -2 * bt;
            out[6] = 0;
            out[7] = 0;
            out[8] = 0;
            out[9] = 0;
            out[10] = 2 * nf;
            out[11] = 0;
            out[12] = (left + right) * lr;
            out[13] = (top + bottom) * bt;
            out[14] = (far + near) * nf;
            out[15] = 1;
        };

        //make the rows into the columns
        static transpose(out, a){
            //If we are transposing ourselves we can skip a few steps but have to cache some values
            if (out === a) {
                var a01 = a[1], a02 = a[2], a03 = a[3], a12 = a[6], a13 = a[7], a23 = a[11];
                out[1] = a[4];
                out[2] = a[8];
                out[3] = a[12];
                out[4] = a01;
                out[6] = a[9];
                out[7] = a[13];
                out[8] = a02;
                out[9] = a12;
                out[11] = a[14];
                out[12] = a03;
                out[13] = a13;
                out[14] = a23;
            }else{
                out[0] = a[0];
                out[1] = a[4];
                out[2] = a[8];
                out[3] = a[12];
                out[4] = a[1];
                out[5] = a[5];
                out[6] = a[9];
                out[7] = a[13];
                out[8] = a[2];
                out[9] = a[6];
                out[10] = a[10];
                out[11] = a[14];
                out[12] = a[3];
                out[13] = a[7];
                out[14] = a[11];
                out[15] = a[15];
            }

            return out;
        }

        //Calculates a 3x3 normal matrix (transpose inverse) from the 4x4 matrix
        static normalMat3(out,a){
            var a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
                a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
                a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
                a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15],

                b00 = a00 * a11 - a01 * a10,
                b01 = a00 * a12 - a02 * a10,
                b02 = a00 * a13 - a03 * a10,
                b03 = a01 * a12 - a02 * a11,
                b04 = a01 * a13 - a03 * a11,
                b05 = a02 * a13 - a03 * a12,
                b06 = a20 * a31 - a21 * a30,
                b07 = a20 * a32 - a22 * a30,
                b08 = a20 * a33 - a23 * a30,
                b09 = a21 * a32 - a22 * a31,
                b10 = a21 * a33 - a23 * a31,
                b11 = a22 * a33 - a23 * a32,

            // Calculate the determinant
            det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

            if (!det) return null;

            det = 1.0 / det;

            out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
            out[1] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
            out[2] = (a10 * b10 - a11 * b08 + a13 * b06) * det;

            out[3] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
            out[4] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
            out[5] = (a01 * b08 - a00 * b10 - a03 * b06) * det;

            out[6] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
            out[7] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
            out[8] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
            return out;
        }

        //New function derived from fromRotationTranslation, just took out the translation stuff.
        static fromQuaternion(out, q){
            // Quaternion math
            var x = q[0], y = q[1], z = q[2], w = q[3],
                x2 = x + x,
                y2 = y + y,
                z2 = z + z,

                xx = x * x2,
                xy = x * y2,
                xz = x * z2,
                yy = y * y2,
                yz = y * z2,
                zz = z * z2,
                wx = w * x2,
                wy = w * y2,
                wz = w * z2;

            out[0] = 1 - (yy + zz);
            out[1] = xy + wz;
            out[2] = xz - wy;
            out[3] = 0;
            out[4] = xy - wz;
            out[5] = 1 - (xx + zz);
            out[6] = yz + wx;
            out[7] = 0;
            out[8] = xz + wy;
            out[9] = yz - wx;
            out[10] = 1 - (xx + yy);
            out[11] = 0;
            return out;
        }

        //https://github.com/toji/gl-matrix/blob/master/src/gl-matrix/mat4.js
        static fromQuaternionTranslation(out, q, v){
            // Quaternion math
            var x = q[0], y = q[1], z = q[2], w = q[3],
                x2 = x + x,
                y2 = y + y,
                z2 = z + z,

                xx = x * x2,
                xy = x * y2,
                xz = x * z2,
                yy = y * y2,
                yz = y * z2,
                zz = z * z2,
                wx = w * x2,
                wy = w * y2,
                wz = w * z2;

            out[0] = 1 - (yy + zz);
            out[1] = xy + wz;
            out[2] = xz - wy;
            out[3] = 0;
            out[4] = xy - wz;
            out[5] = 1 - (xx + zz);
            out[6] = yz + wx;
            out[7] = 0;
            out[8] = xz + wy;
            out[9] = yz - wx;
            out[10] = 1 - (xx + yy);
            out[11] = 0;
            out[12] = v[0];
            out[13] = v[1];
            out[14] = v[2];
            out[15] = 1;
            return out;
        }

        static fromQuaternionTranslationScale(out, q, v, s){
            // Quaternion math
            var x = q[0], y = q[1], z = q[2], w = q[3],
            x2 = x + x,
            y2 = y + y,
            z2 = z + z,

            xx = x * x2,
            xy = x * y2,
            xz = x * z2,
            yy = y * y2,
            yz = y * z2,
            zz = z * z2,
            wx = w * x2,
            wy = w * y2,
            wz = w * z2,
            sx = s[0],
            sy = s[1],
            sz = s[2];

            out[0] = (1 - (yy + zz)) * sx;
            out[1] = (xy + wz) * sx;
            out[2] = (xz - wy) * sx;
            out[3] = 0;
            out[4] = (xy - wz) * sy;
            out[5] = (1 - (xx + zz)) * sy;
            out[6] = (yz + wx) * sy;
            out[7] = 0;
            out[8] = (xz + wy) * sz;
            out[9] = (yz - wx) * sz;
            out[10] = (1 - (xx + yy)) * sz;
            out[11] = 0;
            out[12] = v[0];
            out[13] = v[1];
            out[14] = v[2];
            out[15] = 1;

            return out;
        }

        static getTranslation(out, mat){
            out[0] = mat[12];
            out[1] = mat[13];
            out[2] = mat[14];
            return out;
        }

        static getScaling(out, mat){
            var m11 = mat[0],
                m12 = mat[1],
                m13 = mat[2],
                m21 = mat[4],
                m22 = mat[5],
                m23 = mat[6],
                m31 = mat[8],
                m32 = mat[9],
                m33 = mat[10];
            out[0] = Math.sqrt(m11 * m11 + m12 * m12 + m13 * m13);
            out[1] = Math.sqrt(m21 * m21 + m22 * m22 + m23 * m23);
            out[2] = Math.sqrt(m31 * m31 + m32 * m32 + m33 * m33);
            return out;
        }

        //Returns a quaternion representing the rotational component of a transformation matrix. If a matrix is built with
        //fromRotationTranslation, the returned quaternion will be the same as the quaternion originally supplied
        static getRotation(out, mat){
            // Algorithm taken from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
            var trace = mat[0] + mat[5] + mat[10],
                S = 0;

            if(trace > 0){
                S = Math.sqrt(trace + 1.0) * 2;
                out[3] = 0.25 * S;
                out[0] = (mat[6] - mat[9]) / S;
                out[1] = (mat[8] - mat[2]) / S; 
                out[2] = (mat[1] - mat[4]) / S; 
            }else if( (mat[0] > mat[5]) & (mat[0] > mat[10]) ){ 
                S = Math.sqrt(1.0 + mat[0] - mat[5] - mat[10]) * 2;
                out[3] = (mat[6] - mat[9]) / S;
                out[0] = 0.25 * S;
                out[1] = (mat[1] + mat[4]) / S; 
                out[2] = (mat[8] + mat[2]) / S; 
            }else if(mat[5] > mat[10]){ 
                S = Math.sqrt(1.0 + mat[5] - mat[0] - mat[10]) * 2;
                out[3] = (mat[8] - mat[2]) / S;
                out[0] = (mat[1] + mat[4]) / S; 
                out[1] = 0.25 * S;
                out[2] = (mat[6] + mat[9]) / S; 
            }else{ 
                S = Math.sqrt(1.0 + mat[10] - mat[0] - mat[5]) * 2;
                out[3] = (mat[1] - mat[4]) / S;
                out[0] = (mat[8] + mat[2]) / S;
                out[1] = (mat[6] + mat[9]) / S;
                out[2] = 0.25 * S;
            }
            return out;
        }

        //....................................................................
        //Static Operation

        //https://github.com/gregtatum/mdn-model-view-projection/blob/master/shared/matrices.js
        static multiplyVector(mat4, v) {
            var x = v[0], y = v[1], z = v[2], w = v[3];
            var c1r1 = mat4[ 0], c2r1 = mat4[ 1], c3r1 = mat4[ 2], c4r1 = mat4[ 3],
                c1r2 = mat4[ 4], c2r2 = mat4[ 5], c3r2 = mat4[ 6], c4r2 = mat4[ 7],
                c1r3 = mat4[ 8], c2r3 = mat4[ 9], c3r3 = mat4[10], c4r3 = mat4[11],
                c1r4 = mat4[12], c2r4 = mat4[13], c3r4 = mat4[14], c4r4 = mat4[15];

            return [
                x*c1r1 + y*c1r2 + z*c1r3 + w*c1r4,
                x*c2r1 + y*c2r2 + z*c2r3 + w*c2r4,
                x*c3r1 + y*c3r2 + z*c3r3 + w*c3r4,
                x*c4r1 + y*c4r2 + z*c4r3 + w*c4r4
            ];
        }

        //https://github.com/toji/gl-matrix/blob/master/src/gl-matrix/vec4.js, vec4.transformMat4
        static transformVec4(out, v, m){
            out[0] = m[0] * v[0] + m[4] * v[1] + m[8]	* v[2] + m[12] * v[3];
            out[1] = m[1] * v[0] + m[5] * v[1] + m[9]	* v[2] + m[13] * v[3];
            out[2] = m[2] * v[0] + m[6] * v[1] + m[10]	* v[2] + m[14] * v[3];
            out[3] = m[3] * v[0] + m[7] * v[1] + m[11]	* v[2] + m[15] * v[3];
            return out;
        }

        //From glMatrix
        //Multiple two mat4 together
        static mult(out, a, b){ 
            var a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
                a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
                a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
                a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];

            // Cache only the current line of the second matrix
            var b0  = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
            out[0] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            out[1] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            out[2] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            out[3] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = b[4]; b1 = b[5]; b2 = b[6]; b3 = b[7];
            out[4] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            out[5] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            out[6] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            out[7] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = b[8]; b1 = b[9]; b2 = b[10]; b3 = b[11];
            out[8] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            out[9] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            out[10] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            out[11] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = b[12]; b1 = b[13]; b2 = b[14]; b3 = b[15];
            out[12] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            out[13] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            out[14] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            out[15] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
            return out;	
        }

        //....................................................................
        //Static Transformation
        static scale(out,x,y,z){
            out[0] *= x;
            out[1] *= x;
            out[2] *= x;
            out[3] *= x;
            out[4] *= y;
            out[5] *= y;
            out[6] *= y;
            out[7] *= y;
            out[8] *= z;
            out[9] *= z;
            out[10] *= z;
            out[11] *= z;
            return out;
        };

        static rotateY(out,rad) {
            var s = Math.sin(rad),
                c = Math.cos(rad),
                a00 = out[0],
                a01 = out[1],
                a02 = out[2],
                a03 = out[3],
                a20 = out[8],
                a21 = out[9],
                a22 = out[10],
                a23 = out[11];

            // Perform axis-specific matrix multiplication
            out[0] = a00 * c - a20 * s;
            out[1] = a01 * c - a21 * s;
            out[2] = a02 * c - a22 * s;
            out[3] = a03 * c - a23 * s;
            out[8] = a00 * s + a20 * c;
            out[9] = a01 * s + a21 * c;
            out[10] = a02 * s + a22 * c;
            out[11] = a03 * s + a23 * c;
            return out;
        }

        static rotateX(out,rad) {
            var s = Math.sin(rad),
                c = Math.cos(rad),
                a10 = out[4],
                a11 = out[5],
                a12 = out[6],
                a13 = out[7],
                a20 = out[8],
                a21 = out[9],
                a22 = out[10],
                a23 = out[11];

            // Perform axis-specific matrix multiplication
            out[4] = a10 * c + a20 * s;
            out[5] = a11 * c + a21 * s;
            out[6] = a12 * c + a22 * s;
            out[7] = a13 * c + a23 * s;
            out[8] = a20 * c - a10 * s;
            out[9] = a21 * c - a11 * s;
            out[10] = a22 * c - a12 * s;
            out[11] = a23 * c - a13 * s;
            return out;
        }

        static rotateZ(out,rad){
            var s = Math.sin(rad),
                c = Math.cos(rad),
                a00 = out[0],
                a01 = out[1],
                a02 = out[2],
                a03 = out[3],
                a10 = out[4],
                a11 = out[5],
                a12 = out[6],
                a13 = out[7];

            // Perform axis-specific matrix multiplication
            out[0] = a00 * c + a10 * s;
            out[1] = a01 * c + a11 * s;
            out[2] = a02 * c + a12 * s;
            out[3] = a03 * c + a13 * s;
            out[4] = a10 * c - a00 * s;
            out[5] = a11 * c - a01 * s;
            out[6] = a12 * c - a02 * s;
            out[7] = a13 * c - a03 * s;
            return out;
        }

        static rotate(out, rad, axis){
            var x = axis[0], y = axis[1], z = axis[2],
                len = Math.sqrt(x * x + y * y + z * z),
                s, c, t,
                a00, a01, a02, a03,
                a10, a11, a12, a13,
                a20, a21, a22, a23,
                b00, b01, b02,
                b10, b11, b12,
                b20, b21, b22;

            if (Math.abs(len) < 0.000001) { return null; }

            len = 1 / len;
            x *= len;
            y *= len;
            z *= len;

            s = Math.sin(rad);
            c = Math.cos(rad);
            t = 1 - c;

            a00 = out[0]; a01 = out[1]; a02 = out[2]; a03 = out[3];
            a10 = out[4]; a11 = out[5]; a12 = out[6]; a13 = out[7];
            a20 = out[8]; a21 = out[9]; a22 = out[10]; a23 = out[11];

            // Construct the elements of the rotation matrix
            b00 = x * x * t + c; b01 = y * x * t + z * s; b02 = z * x * t - y * s;
            b10 = x * y * t - z * s; b11 = y * y * t + c; b12 = z * y * t + x * s;
            b20 = x * z * t + y * s; b21 = y * z * t - x * s; b22 = z * z * t + c;

            // Perform rotation-specific matrix multiplication
            out[0] = a00 * b00 + a10 * b01 + a20 * b02;
            out[1] = a01 * b00 + a11 * b01 + a21 * b02;
            out[2] = a02 * b00 + a12 * b01 + a22 * b02;
            out[3] = a03 * b00 + a13 * b01 + a23 * b02;
            out[4] = a00 * b10 + a10 * b11 + a20 * b12;
            out[5] = a01 * b10 + a11 * b11 + a21 * b12;
            out[6] = a02 * b10 + a12 * b11 + a22 * b12;
            out[7] = a03 * b10 + a13 * b11 + a23 * b12;
            out[8] = a00 * b20 + a10 * b21 + a20 * b22;
            out[9] = a01 * b20 + a11 * b21 + a21 * b22;
            out[10] = a02 * b20 + a12 * b21 + a22 * b22;
            out[11] = a03 * b20 + a13 * b21 + a23 * b22;
        }

        static invert(out,mat) {
            if(mat === undefined) mat = out; //If input isn't sent, then output is also input

            var a00 = mat[0], a01 = mat[1], a02 = mat[2], a03 = mat[3],
                a10 = mat[4], a11 = mat[5], a12 = mat[6], a13 = mat[7],
                a20 = mat[8], a21 = mat[9], a22 = mat[10], a23 = mat[11],
                a30 = mat[12], a31 = mat[13], a32 = mat[14], a33 = mat[15],

                b00 = a00 * a11 - a01 * a10,
                b01 = a00 * a12 - a02 * a10,
                b02 = a00 * a13 - a03 * a10,
                b03 = a01 * a12 - a02 * a11,
                b04 = a01 * a13 - a03 * a11,
                b05 = a02 * a13 - a03 * a12,
                b06 = a20 * a31 - a21 * a30,
                b07 = a20 * a32 - a22 * a30,
                b08 = a20 * a33 - a23 * a30,
                b09 = a21 * a32 - a22 * a31,
                b10 = a21 * a33 - a23 * a31,
                b11 = a22 * a33 - a23 * a32,

                // Calculate the determinant
                det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

            if (!det) return false;
            det = 1.0 / det;

            out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
            out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
            out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
            out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
            out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
            out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
            out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
            out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
            out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
            out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
            out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
            out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
            out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
            out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
            out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
            out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;

            return true;
        }

        //https://github.com/toji/gl-matrix/blob/master/src/gl-matrix/mat4.js  mat4.scalar.translate = function (out, a, v) {
        static translate(out,x,y,z){
            out[12] = out[0] * x + out[4] * y + out[8]	* z + out[12];
            out[13] = out[1] * x + out[5] * y + out[9]	* z + out[13];
            out[14] = out[2] * x + out[6] * y + out[10]	* z + out[14];
            out[15] = out[3] * x + out[7] * y + out[11]	* z + out[15];
        }
    //endregion
}

function mix(a, b, t) {
    return [(1 - t) * a[0] + t * b[0], (1 - t) * a[1] + t * b[1]];
}

// Fungsi Bezier custom
function cubicMix(a, b, c, d, t) {
    var ab = mix(a, b, t);
    var bc = mix(b, c, t);
    var cd = mix(c, d, t);
    var abc = mix(ab, bc, t);
    var bcd = mix(bc, cd, t);
    return mix(abc, bcd, t);
}

function normalizeColor(rgb) {
    if (rgb.length !== 3) {
        throw new Error("Input must be a list of 3 integers representing RGB values (0-255).");
    }
    return rgb.map(channel => channel / 255.0);
}

function translate(m, mTrans, normal = [0,0,0]){
    for(var i = 0; i < m.length; i+=3){
        m[i] += mTrans[0]+normal[0];
        m[i+1] += mTrans[1]+normal[1];
        m[i+2] += mTrans[2]+normal[2];
    }
}

function scale(m, scale){
    for (var i = 0; i < m.length; i += 3) {
        m[i] *= scale[0] > 0 ? scale[0] : 1; // Apply scale only if positive
        m[i + 1] *= scale[1] > 0 ? scale[1] : 1;
        m[i + 2] *= scale[2] > 0 ? scale[2] : 1;
    }
}

function rotate(points, yDeg, zDeg, xDeg) {
    var yDeg = Math.PI * yDeg / 180; // Convert angle from degrees to radians
    var zDeg = Math.PI * zDeg / 180; // Convert angle from degrees to radians
    var xDeg = Math.PI * xDeg / 180; // Convert angle from degrees to radians
    var cosa = Math.cos(xDeg);
    var sina = Math.sin(xDeg);

    var cosb = Math.cos(yDeg);
    var sinb = Math.sin(yDeg);

    var cosc = Math.cos(zDeg);
    var sinc = Math.sin(zDeg);

    var Axx = cosa*cosb;
    var Axy = cosa*sinb*sinc - sina*cosc;
    var Axz = cosa*sinb*cosc + sina*sinc;

    var Ayx = sina*cosb;
    var Ayy = sina*sinb*sinc + cosa*cosc;
    var Ayz = sina*sinb*cosc - cosa*sinc;

    var Azx = -sinb;
    var Azy = cosb*sinc;
    var Azz = cosb*cosc;

    for (var i = 0; i < points.length; i+= 3) {
        var px = points[i];
        var py = points[i+1];
        var pz = points[i+2];

        points[i] = Axx*px + Axy*py + Axz*pz;
        points[i+1] = Ayx*px + Ayy*py + Ayz*pz;
        points[i+2] = Azx*px + Azy*py + Azz*pz;
    }
}

function generateBSpline(controlPoint, m, degree){
    var curves = [];
    var knotVector = []

    var n = controlPoint.length/2;


    // Calculate the knot values based on the degree and number of control points
    for (var i = 0; i < n + degree+1; i++) {
        if (i < degree + 1) {
        knotVector.push(0);
        } else if (i >= n) {
        knotVector.push(n - degree);
        } else {
        knotVector.push(i - degree);
        }
    }



    var basisFunc = function(i,j,t){
        if (j == 0){
            if(knotVector[i] <= t && t<(knotVector[(i+1)])){ 
            return 1;
            }else{
            return 0;
            }
        }

        var den1 = knotVector[i + j] - knotVector[i];
        var den2 = knotVector[i + j + 1] - knotVector[i + 1];
        
        var term1 = 0;
        var term2 = 0;
        

        if (den1 != 0 && !isNaN(den1)) {
            term1 = ((t - knotVector[i]) / den1) * basisFunc(i,j-1,t);
        }
        
        if (den2 != 0 && !isNaN(den2)) {
            term2 = ((knotVector[i + j + 1] - t) / den2) * basisFunc(i+1,j-1,t);
        }
        
        return term1 + term2;
    }


    for(var t=0;t<m;t++){
        var x=0;
        var y=0;
        
        var u = (t/m * (knotVector[controlPoint.length/2] - knotVector[degree]) ) + knotVector[degree] ;

        //C(t)
        for(var key =0;key<n;key++){

        var C = basisFunc(key,degree,u);
        console.log(C);
        x+=(controlPoint[key*2] * C);
        y+=(controlPoint[key*2+1] * C);
        console.log(t+" "+degree+" "+x+" "+y+" "+C);
        }
        curves.push(x);
        curves.push(y);
        
    }
    console.log(curves)
    return curves;
}

function generateColors(vertices, color, isRainbow = false) {
    var colors = [];
    var definedColors = [[.1, .0, .0, 1],
    [.0, .1, .0, 1],
    [.0, .0, .1, 1]]

    if(isRainbow){
        for (var i = 0; i < vertices.length; i += 4) {
            colors[i] = definedColors[colors.length % definedColors.length];
        }
        return colors.reduce(function (a, b) { return a.concat(b); });
    }

    for (var i = 0; i < vertices.length; i += 3) {
        colors.push(color[0], color[1], color[2]);
    }
    return colors;
};

function main(type = 1){
    var CANVAS = document.getElementById("mycanvas");
    
    CANVAS.width = window.innerWidth;
    CANVAS.height = window.innerHeight;

    // EVENT LISTENER
    var AMORTIZATION = 0.95;
    var x_prev, y_prev;
    var drag = false;
    var THETA = 0, PHI = 0;
    var dX = 0, dY = 0;
    var zoomFactor = 1;

    var mouseDown = function(e){
        drag = true;
        x_prev = e.pageX;
        y_prev = e.pageY;

        // BIAR GABISA DIKLIK KANAN
        // e.preventDefault();

        return false;
    };
    
    var mouseUp = function(e){
        drag = false;
    };

    var mouseMove = function(e){
        if(!drag) return false;

        dX = (e.pageX - x_prev) * 2 * Math.PI / CANVAS.width;
        dY = (e.pageY - y_prev) * 2 * Math.PI / CANVAS.height;

        THETA += dX;
        PHI += dY;

        x_prev = e.pageX;
        y_prev = e.pageY;

        // BIAR GABISA DIKLIK KANAN
        // e.preventDefault();

        return false;
    };

    var keyDown = function(e) {
        var speed = 0.01; // Kecepatan pergerakan, ubah sesuai kebutuhan
        var zoomSpeed = 0.1;
    
        switch(e.key.toLowerCase()) {
            case 'w': // Keatas
                PHI -= speed;
                break;
            case 'a': // Kekiri
                THETA -= speed;
                break;
            case 's': // Kebawah
                PHI += speed;
                break;
            case 'd': // Kekanan
                THETA += speed;
                break;
            case 'q': // Zoom out
                zoomFactor -= zoomSpeed;
                break;
            case 'e': // Zoom in
                zoomFactor += zoomSpeed;
                break;
        }
    
        // BIAR GABISA DIKLIK KANAN
        // e.preventDefault();
    
        return false;
    }
    
    CANVAS.addEventListener("mousedown", mouseDown, false);
    CANVAS.addEventListener("mouseup", mouseUp, false);
    CANVAS.addEventListener("mouseout", mouseUp, false);
    CANVAS.addEventListener("mousemove", mouseMove, false);
    
    window.addEventListener("keydown", keyDown, false);

    var GL;
    try{
        GL = CANVAS.getContext("webgl", {antialias: true});
    }catch(e){
        alert("WebGL context cannot be initialized");
        return false;
    }


    //shaders
    // Vertex Shader
//     var shader_vertex_source2d = `
//     attribute vec3 a_position0; // Control point 0
//     attribute vec3 a_position1; // Control point 1
//     attribute vec3 a_position2; // Control point 2
//     attribute vec3 a_position3; // Control point 3

//     uniform mat4 PMatrix;
//     uniform mat4 VMatrix;
//     uniform mat4 MMatrix;
    
//     uniform float u_t; // Parameter (0.0 to 1.0)
    
//     void main() {
//       // Calculate the Bezier curve position using cubic interpolation
//       vec3 position = (1.0 - u_t) * (1.0 - u_t) * (1.0 - u_t) * a_position0 +
//                      3.0 * u_t * (1.0 - u_t) * (1.0 - u_t) * a_position1 +
//                      3.0 * u_t * u_t * (1.0 - u_t) * a_position2 +
//                      u_t * u_t * u_t * a_position3;
    
//       // Transform the position to clip space
//       gl_Position = PMatrix * VMatrix * MMatrix * vec4(position, 1.0);
//     }

//   `;
    var shader_vertex_source=`
    attribute vec3 position;
    attribute vec3 color;

    uniform mat4 PMatrix;
    uniform mat4 VMatrix;
    uniform mat4 MMatrix;
    
    varying vec3 vColor;
    void main(void) {
    gl_Position = PMatrix*VMatrix*MMatrix*vec4(position, 1.);
    gl_PointSize = 20.0;
    vColor = color;
    }`;
    var shader_fragment_source =`
    precision mediump float;
    varying vec3 vColor;
    // uniform vec3 color;
    void main(void) {
    gl_FragColor = vec4(vColor, 1.);
    
    }`;
    var compile_shader = function(source, type, typeString) {
        var shader = GL.createShader(type);
        GL.shaderSource(shader, source);
        GL.compileShader(shader);
        if (!GL.getShaderParameter(shader, GL.COMPILE_STATUS)) {
        alert("ERROR IN " + typeString + " SHADER: " + GL.getShaderInfoLog(shader));
        return false;
        }
        return shader;
    };
    
    // 2D
    // var shader_vertex2D = compile_shader(shader_vertex_source2d, GL.VERTEX_SHADER, "VERTEX");
    
    var shader_fragment = compile_shader(shader_fragment_source, GL.FRAGMENT_SHADER, "FRAGMENT");
    
    // var SHADER_PROGRAM2D = GL.createProgram();
    // GL.attachShader(SHADER_PROGRAM2D, shader_vertex2D);
    
    // GL.linkProgram(SHADER_PROGRAM2D);
    
    var shader_vertex = compile_shader(shader_vertex_source, GL.VERTEX_SHADER, "VERTEX");
    
    var SHADER_PROGRAM = GL.createProgram();
    GL.attachShader(SHADER_PROGRAM, shader_vertex);
    GL.attachShader(SHADER_PROGRAM, shader_fragment);
    
    GL.linkProgram(SHADER_PROGRAM);


    var _color = GL.getAttribLocation(SHADER_PROGRAM, "color");
    var _position = GL.getAttribLocation(SHADER_PROGRAM, "position");
    var u_colorLocation = GL.getUniformLocation(SHADER_PROGRAM, "u_color");
    var position_vao = GL.getAttribLocation(SHADER_PROGRAM, "a_position0");
    var controlPoint1_vao = GL.getAttribLocation(SHADER_PROGRAM, "a_position1");
    var controlPoint2_vao = GL.getAttribLocation(SHADER_PROGRAM, "a_position2");
    var endPoint_vao = GL.getAttribLocation(SHADER_PROGRAM, "a_position3");
    var u_tLocation = GL.getUniformLocation(SHADER_PROGRAM, "u_t");
    

    
    var _PMatrix = GL.getUniformLocation(SHADER_PROGRAM,"PMatrix"); //projection
    var _VMatrix = GL.getUniformLocation(SHADER_PROGRAM,"VMatrix"); //View
    var _MMatrix = GL.getUniformLocation(SHADER_PROGRAM,"MMatrix"); //Model


    GL.enableVertexAttribArray(_color);
    GL.enableVertexAttribArray(_position);
    GL.enableVertexAttribArray(position_vao);
    GL.enableVertexAttribArray(controlPoint1_vao);
    GL.enableVertexAttribArray(controlPoint2_vao);
    GL.enableVertexAttribArray(endPoint_vao);
    GL.useProgram(SHADER_PROGRAM);

    var arrColors = [
        [1, 0, 0], // merah
        [0.3, .8, 0.5], // hijau
        [0.5, .8, 0.5], // DARKER GREEN
        [1, 1, 0], // kuning
        [1, 0, 1], // magenta
        [0, 1, 1], // cyan
        [1, 1, 1], // putih
        [0.42, 0.286, 0.169], // brown
        [0, 0, 0],  // hitam
        [255,218,161], //kuning pooh agak gelap
        [255,241,156], //kuning pooh
        [255,109,109], //merah pooh
        [216,255,254], //biru pooh
        [255,255,224], //putih pooh
        [255, 255,103], //Main Rabbit
        [255, 153,156], // Second Rabbit
    ];

    var arrColors2 = [
        [1, 0, 0], // merah
        [0.3, .8, 0.5], // hijau
        [0.5, .8, 0.5], // DARKER GREEN
        [1, 1, 0], // kuning
        [1, 0, 1], // magenta
        [0, 1, 1], // cyan
        [1, 1, 1], // putih
        [0.42, 0.286, 0.169], // brown
        [0, 0, 0],  // hitam
        [255,218,161], //kuning pooh agak gelap
        [255,241,156], //kuning pooh
        [255,109,109], //merah pooh
        [216,255,254], //biru pooh
        [255,255,224], //putih pooh
        [255,204,204], //pink cerah piglet
        [219, 112, 147], //pink piglet
        [253,82,167], //pink gelap badan piglet
        [255, 255, 102], //matahari bulat
        [102, 204,0] //surface
    ];

    var controlPoints = [
        // Top left corner (A)
        [5, 12],
    
        // Top control points (additional)
        [6, 11.5],
        [7, 11],
    
        // Top middle control point (A1) - Adjust y-coordinate for more curvature
        [8, 10],
    
        // Top control points (additional)
        [9, 11],
        [10, 11.5],
    
        // Top right corner (B)
        [15, 12],
    
        // Bottom right corner (C)
        [15, 7],
    
        // Bottom control points (additional)
        [14, 7.5],
        [13, 8],
    
        // Bottom middle control point (C1) - Adjust y-coordinate for more curvature
        [8, 6],
    
        // Bottom control points (additional)
        [7, 8],
        [6, 7.5],
    
        // Bottom left corner (D)
        [5, 7],
    ];

    // POSITION POOH
    var vecTransPooh = [10,0,-1]
    
    // POSITION Rabbit
    var vecTransRabbit = [-12,0,-1]

    var translationVectorPiglet = [5, 0, 1]; // [dx, dy, dz] for 3D translation
    var translationVectorPohon = [-6,-5.5,-8.3]
    var translationVectorSun = [0.5,-4,-1]; //kanan kiri,depan blakang, atas bawah
    var translationVectorCloud = [14,-5,1];
    var translationVectorCloud2 = [-16,5,1];

   

    //matrix
    var PROJECTION_MATRIX = LIBS.get_projection(40, CANVAS.width/CANVAS.height, 1,100);
    var VIEW_MATRIX = LIBS.get_I4();
    var MODEL_MATRIX = LIBS.get_I4();

    LIBS.translateZ(VIEW_MATRIX,-80);


    /*========================= DRAWING ========================= */
    GL.clearColor(0.0, 0.0, 0.0, 0.0);


    GL.enable(GL.DEPTH_TEST);
    GL.depthFunc(GL.LEQUAL);
    

    var count,y,z;
    var time_prev = 0;
    var animate = function(time) {
        if(time > 0){
            var dt = time - time_prev;
            
            if(!drag){
                dX *= AMORTIZATION;
                dY *= AMORTIZATION;
                THETA += dX;
                PHI += dY;
            }
            
            LIBS.set_I4(MODEL_MATRIX);
            LIBS.scale(MODEL_MATRIX, [zoomFactor, zoomFactor, zoomFactor])
            // LIBS.rotateZ(MODEL_MATRIX, dt * LIBS.degToRad(0.1));
            LIBS.rotateY(MODEL_MATRIX, THETA + LIBS.degToRad(0));
            LIBS.rotateX(MODEL_MATRIX, PHI + LIBS.degToRad(110));
            // LIBS.rotateY(MODEL_MATRIX, dt * LIBS.degToRad(0.1));
            // LIBS.rotateX(MODEL_MATRIX, dt * LIBS.degToRad(0.1));

            // ANIMASI POOH
            y = dt*.001;
            z = -dt*.003
            translate(vecTransRabbit, [0, y, z])
            translate(vecTransPooh, [0, y, z])
            translate(translationVectorPiglet, [0, y, z])


            time_prev = time;
        }
        GL.viewport(0, 0, CANVAS.width, CANVAS.height);
        GL.clear(GL.COLOR_BUFFER_BIT | GL.DEPTH_BUFFER_BIT);

         // SUN
    var sun = new Icosahedron3D(1);
    sun_vertices = sun.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    sun_indices = sun.TriangleIndices;
    sun_color = generateColors(sun_vertices, normalizeColor(arrColors2[17]), false);

    scale(sun_vertices, [2,2,2])
    // rotate(head_vertices, 0, 0, 90)
    translate(sun_vertices, [-8+translationVectorSun[0],-8+translationVectorSun[1],-15+translationVectorSun[2]])

    var vertex_bufferSUN = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferSUN);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(sun_vertices), GL.STATIC_DRAW);

    var color_bufferSUN = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferSUN);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(sun_color), GL.STATIC_DRAW);

    var index_bufferSUN = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferSUN);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(sun_indices), GL.STATIC_DRAW);

    // CLOUD BULAT 1
    var cloudB1 = new Icosahedron3D(1);
    cloudB1_vertices = cloudB1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    cloudB1_indices = cloudB1.TriangleIndices;
    cloudB1_color = generateColors(cloudB1_vertices, normalizeColor(arrColors2[12]), false);

    scale(cloudB1_vertices, [4,2,2])
    // rotate(head_vertices, 0, 0, 90)
    translate(cloudB1_vertices, [-8+translationVectorCloud[0],-8+translationVectorCloud[1],-15+translationVectorCloud[2]])

    var vertex_buffercloudB1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB1_vertices), GL.STATIC_DRAW);

    var color_buffercloudB1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB1_color), GL.STATIC_DRAW);

    var index_buffercloudB1 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB1);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB1_indices), GL.STATIC_DRAW);

    // CLOUD BULAT 2
    var cloudB2 = new Icosahedron3D(1);
    cloudB2_vertices = cloudB2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    cloudB2_indices = cloudB2.TriangleIndices;
    cloudB2_color = generateColors(cloudB2_vertices, normalizeColor(arrColors2[12]), false);

    scale(cloudB2_vertices, [4,2.5,2.5])
    // rotate(head_vertices, 0, 0, 90)
    translate(cloudB2_vertices, [-2+translationVectorCloud[0],-8+translationVectorCloud[1],-15+translationVectorCloud[2]])

    var vertex_buffercloudB2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB2_vertices), GL.STATIC_DRAW);

    var color_buffercloudB2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB2_color), GL.STATIC_DRAW);

    var index_buffercloudB2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB2_indices), GL.STATIC_DRAW);

     // CLOUD BULAT 3
     var cloudB3 = new Icosahedron3D(1);
     cloudB3_vertices = cloudB3.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
     cloudB3_indices = cloudB3.TriangleIndices;
     cloudB3_color = generateColors(cloudB3_vertices, normalizeColor(arrColors2[12]), false);
 
     scale(cloudB3_vertices, [4,2.5,2.5])
     // rotate(head_vertices, 0, 0, 90)
     translate(cloudB3_vertices, [-4+translationVectorCloud[0],-8+translationVectorCloud[1],-16+translationVectorCloud[2]])
 
     var vertex_buffercloudB3 = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB3);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB3_vertices), GL.STATIC_DRAW);
 
     var color_buffercloudB3 = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB3);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB3_color), GL.STATIC_DRAW);
 
     var index_buffercloudB3 = GL.createBuffer();
     GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB3);
     GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB3_indices), GL.STATIC_DRAW);

    
    // CLOUD BULAT 4
    var cloudB4 = new Icosahedron3D(1);
    cloudB4_vertices = cloudB4.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    cloudB4_indices = cloudB4.TriangleIndices;
    cloudB4_color = generateColors(cloudB4_vertices, normalizeColor(arrColors2[12]), false);

    scale(cloudB4_vertices, [4,2,2])
    // rotate(head_vertices, 0, 0, 90)
    translate(cloudB4_vertices, [-8+translationVectorCloud2[0],-8+translationVectorCloud2[1],-15+translationVectorCloud2[2]])

    var vertex_buffercloudB4 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB4);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB4_vertices), GL.STATIC_DRAW);

    var color_buffercloudB4 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB4);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB4_color), GL.STATIC_DRAW);

    var index_buffercloudB4 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB4);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB4_indices), GL.STATIC_DRAW);

    // CLOUD BULAT 5
    var cloudB5 = new Icosahedron3D(1);
    cloudB5_vertices = cloudB5.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    cloudB5_indices = cloudB5.TriangleIndices;
    cloudB5_color = generateColors(cloudB5_vertices, normalizeColor(arrColors2[12]), false);

    scale(cloudB5_vertices, [4,2.5,2.5])
    // rotate(head_vertices, 0, 0, 90)
    translate(cloudB5_vertices, [-2+translationVectorCloud2[0],-8+translationVectorCloud2[1],-15+translationVectorCloud2[2]])

    var vertex_buffercloudB5 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB5);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB5_vertices), GL.STATIC_DRAW);

    var color_buffercloudB5 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB5);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB5_color), GL.STATIC_DRAW);

    var index_buffercloudB5 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB5);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB5_indices), GL.STATIC_DRAW);

     // CLOUD BULAT 6
     var cloudB6 = new Icosahedron3D(1);
     cloudB6_vertices = cloudB6.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
     cloudB6_indices = cloudB6.TriangleIndices;
     cloudB6_color = generateColors(cloudB6_vertices, normalizeColor(arrColors2[12]), false);
 
     scale(cloudB6_vertices, [4,2.5,2.5])
     // rotate(head_vertices, 0, 0, 90)
     translate(cloudB6_vertices, [-4+translationVectorCloud2[0],-8+translationVectorCloud2[1],-16+translationVectorCloud2[2]])
 
     var vertex_buffercloudB6 = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB6);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB6_vertices), GL.STATIC_DRAW);
 
     var color_buffercloudB6 = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB6);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cloudB6_color), GL.STATIC_DRAW);
 
     var index_buffercloudB6 = GL.createBuffer();
     GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB6);
     GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(cloudB6_indices), GL.STATIC_DRAW);

    
    // SURFACE
    var objSurface = generateCube(0, 0, 4, [35,35,1], normalizeColor(arrColors2[18]), false);
    var SURFACE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, SURFACE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objSurface.vertices), GL.STATIC_DRAW);
    
    var SURFACE_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, SURFACE_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objSurface.colors), GL.STATIC_DRAW);

    var SURFACE_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, SURFACE_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(objSurface.faces), GL.STATIC_DRAW);
    
    // TBRANCH
    var objTb = generateCylinder(0, 0, -0.25, 4, 1000, arrColors2[7], false);
    scale(objTb.vertices, [1.2,1,-5])
    translate(objTb.vertices, [0+translationVectorPohon[0],0+translationVectorPohon[1],0+translationVectorPohon[2]])
    var TBRANCH2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objTb.vertices), GL.STATIC_DRAW);
    
    var TBRANCH2_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH2_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objTb.colors), GL.STATIC_DRAW);

    var TBRANCH2_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TBRANCH2_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(objTb.faces), GL.STATIC_DRAW);
    
    // TBRANCH
    var objTb2 = generateCylinder(0, 0, -0.25, 7.2, 1000, arrColors2[7], false);
    scale(objTb2.vertices, [1.2,1,-1])
    translate(objTb2.vertices, [0+translationVectorPohon[0],1+translationVectorPohon[1],5+translationVectorPohon[2]])
    var TBRANCH_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objTb2.vertices), GL.STATIC_DRAW);
    
    var TBRANCH_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objTb2.colors), GL.STATIC_DRAW);

    var TBRANCH_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TBRANCH_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(objTb2.faces), GL.STATIC_DRAW);
    
    // TLEAF
    var objLeaf = generateCone(0, 0, -9, 2.4, 0.6, 100, 
                                arrColors2[2], false);
    translate(objLeaf.vertices, [0+translationVectorPohon[0],1+translationVectorPohon[1],5+translationVectorPohon[2]])
    var TLEAF_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objLeaf.vertices), GL.STATIC_DRAW);
    
    var TLEAF_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objLeaf.colors), GL.STATIC_DRAW);

    var TLEAF_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TLEAF_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(objLeaf.faces), GL.STATIC_DRAW);


    // TLEAF
    var objLeaf2 = generateCone2(0, 1, 0, 2.4, 0.3, 100, 
                                arrColors2[1], false);
    translate(objLeaf2.vertices, [0+translationVectorPohon[0],1+translationVectorPohon[1],5+translationVectorPohon[2]])
    var TLEAF2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objLeaf2.vertices), GL.STATIC_DRAW);

    var TLEAF2_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF2_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(objLeaf2.colors), GL.STATIC_DRAW);

    var TLEAF2_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TLEAF2_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(objLeaf2.faces), GL.STATIC_DRAW);





    // HEAD PIGLET
    var icosahedronP1 = new Icosahedron3D(1);
    head_verticesP1 = icosahedronP1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    head_indicesP1 = icosahedronP1.TriangleIndices;
    head_colorP1 = generateColors(head_verticesP1, normalizeColor(arrColors2[14]), false);

     scale(head_verticesP1, [0,1,1.5])
    translate(head_verticesP1, [0+translationVectorPiglet[0],0+translationVectorPiglet[1],-3+translationVectorPiglet[2]])
    // rotate(head_vertices, 0, 0, 90)

    var vertex_bufferP1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_verticesP1), GL.STATIC_DRAW);

    var color_bufferP1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_colorP1), GL.STATIC_DRAW);

    var index_bufferP1 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP1);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(head_indicesP1), GL.STATIC_DRAW);

    // BODY PIGLET
    var icosahedronP2 = new Icosahedron3D(1);
    head_verticesP2 = icosahedronP2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    head_indicesP2 = icosahedronP2.TriangleIndices;
    head_colorP2 = generateColors(head_verticesP2, normalizeColor(arrColors2[16]), false);

    scale(head_verticesP2, [1.2,1,2])
    translate(head_verticesP2, [0+translationVectorPiglet[0],0+translationVectorPiglet[1],0+translationVectorPiglet[2]])
    // rotate(head_vertices, 0, 0, 90)

    var vertex_bufferP2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_verticesP2), GL.STATIC_DRAW);

    var color_bufferP2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_colorP2), GL.STATIC_DRAW);

    var index_bufferP2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(head_indicesP2), GL.STATIC_DRAW);

    // NOSE PIGLET
    var nose = generateCylinder(0, 0, -0.25, 4, 1000, normalizeColor(arrColors2[15]), false);
    scale(nose.vertices, [1,-14,0.2])
    rotate(nose.vertices, 0, 90, 180)
    translate(nose.vertices, [0+translationVectorPiglet[0],.35+translationVectorPiglet[1],-2.6+translationVectorPiglet[2]])
    var NOSE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose.vertices), GL.STATIC_DRAW);
    
    var NOSE_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose.colors), GL.STATIC_DRAW);

    var NOSE_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(nose.faces), GL.STATIC_DRAW);

    // NOSE PIGLET BULAT
    var icosahedronP5 = new Icosahedron3D(1);
    rnose_verticesP1 = icosahedronP5.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rnose_indicesP1 = icosahedronP5.TriangleIndices;
    rnose_colorP1 = generateColors(rnose_verticesP1, normalizeColor(arrColors2[15]), false);

    scale(rnose_verticesP1, [0.26,0.26,0.26])
    translate(rnose_verticesP1, [0+translationVectorPiglet[0],1.11+translationVectorPiglet[1],-2.6+translationVectorPiglet[2]])
    // rotate(head_vertices, 0, 0, 90)

    var vertex_bufferP5 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP5);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rnose_verticesP1), GL.STATIC_DRAW);

    var color_bufferP5 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP5);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rnose_colorP1), GL.STATIC_DRAW);

    var index_bufferP5 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP5);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rnose_indicesP1), GL.STATIC_DRAW);


    // MATA KIRI
    var icosahedronP3 = new Icosahedron3D(1);
    lefteye_verticesP1 = icosahedronP3.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    lefteye_indicesP1 = icosahedronP3.TriangleIndices;
    lefteye_colorP1 = generateColors(lefteye_verticesP1, normalizeColor(arrColors2[8]), false);

    scale(lefteye_verticesP1, [0.09,0.09,0.09])
    translate(lefteye_verticesP1, [-0.35+translationVectorPiglet[0],0.95+translationVectorPiglet[1],-3.25+translationVectorPiglet[2]])
    // rotate(head_vertices, 0, 0, 90)

    var vertex_bufferP3 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP3);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefteye_verticesP1), GL.STATIC_DRAW);

    var color_bufferP3 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP3);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefteye_colorP1), GL.STATIC_DRAW);

    var index_bufferP3 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP3);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(lefteye_indicesP1), GL.STATIC_DRAW);

    // MATA KANAN
    var icosahedronP4 = new Icosahedron3D(1);
    righteye_verticesP1 = icosahedronP4.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    righteye_indicesP1 = icosahedronP4.TriangleIndices;
    righteye_colorP1 = generateColors(righteye_verticesP1, normalizeColor(arrColors2[8]), false);

    scale(righteye_verticesP1, [0.09,0.09,0.09])
    translate(righteye_verticesP1, [0.35+translationVectorPiglet[0],0.95+translationVectorPiglet[1],-3.25+translationVectorPiglet[2]])
    // rotate(head_vertices, 0, 0, 90)

    var vertex_bufferP4 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP4);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righteye_verticesP1), GL.STATIC_DRAW);

    var color_bufferP4 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP4);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righteye_colorP1), GL.STATIC_DRAW);

    var index_bufferP4 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP4);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(righteye_indicesP1), GL.STATIC_DRAW);

    // TELINGA KIRI
    var L_ear = generateCone(0, 0, -0.25, 0.5, 5, 1000, normalizeColor(arrColors2[15]), false);
    scale(L_ear.vertices, [1.2,0.5,0.30])
    rotate(L_ear.vertices, -45, 180, 180)
    translate(L_ear.vertices, [-0.5+translationVectorPiglet[0],.22+translationVectorPiglet[1],-4+translationVectorPiglet[2]])
    var L_EAR_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, L_EAR_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(L_ear.vertices), GL.STATIC_DRAW);
    
    var L_EAR_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, L_EAR_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(L_ear.colors), GL.STATIC_DRAW);

    var L_EAR_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, L_EAR_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(L_ear.faces), GL.STATIC_DRAW);
    
    // TELINGA KANAN
    var R_ear = generateCone(0, 0, -0.25, 0.5, 5, 1000, normalizeColor(arrColors2[15]), false);
    scale(R_ear.vertices, [1.2,0.5,0.30])
    rotate(R_ear.vertices, 45, 180, 180)
    translate(R_ear.vertices, [0.5+translationVectorPiglet[0],.22+translationVectorPiglet[1],-4+translationVectorPiglet[2]])
    var R_EAR_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, R_EAR_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(R_ear.vertices), GL.STATIC_DRAW);
    
    var R_EAR_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, R_EAR_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(R_ear.colors), GL.STATIC_DRAW);

    var R_EAR_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, R_EAR_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(R_ear.faces), GL.STATIC_DRAW);

    // TANGAN KIRI ATAS
    var icosahedronP6 = new Icosahedron3D(1);
    upperleft_verticesP1 = icosahedronP6.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    upperleft_indicesP1 = icosahedronP6.TriangleIndices;
    upperleft_colorP1 = generateColors(upperleft_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(upperleft_verticesP1, [0.5,0.35,1.1])
    rotate(upperleft_verticesP1, 0, -40, 90)
    translate(upperleft_verticesP1, [-1+translationVectorPiglet[0],0+translationVectorPiglet[1],-0.85+translationVectorPiglet[2]])

    var vertex_bufferP6 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP6);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(upperleft_verticesP1), GL.STATIC_DRAW);

    var color_bufferP6 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP6);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(upperleft_colorP1), GL.STATIC_DRAW);

    var index_bufferP6 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP6);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(upperleft_indicesP1), GL.STATIC_DRAW);

    // TANGAN KIRI BAWAH
    var icosahedronP7 = new Icosahedron3D(1);
    bottomleft_verticesP1 = icosahedronP7.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    bottomleft_indicesP1 = icosahedronP7.TriangleIndices;
    bottomleft_colorP1 = generateColors(bottomleft_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(bottomleft_verticesP1, [0.5,0.35,1.1])
    //x y z -> z = naik turun, y = kiri kanan, x = depan belakang
    rotate(bottomleft_verticesP1, 0, 40, 90)
    translate(bottomleft_verticesP1, [-2.15+translationVectorPiglet[0],0+translationVectorPiglet[1],-0.85+translationVectorPiglet[2]])

    var vertex_bufferP7 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP7);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bottomleft_verticesP1), GL.STATIC_DRAW);

    var color_bufferP7 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP7);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bottomleft_colorP1), GL.STATIC_DRAW);

    var index_bufferP7 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP7);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(bottomleft_indicesP1), GL.STATIC_DRAW);

    // TANGAN KANAN ATAS
    var icosahedronP8 = new Icosahedron3D(1);
    upperright_verticesP1 = icosahedronP8.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    upperright_indicesP1 = icosahedronP8.TriangleIndices;
    upperright_colorP1 = generateColors(upperright_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(upperright_verticesP1, [0.5,0.35,1.1])
    rotate(upperright_verticesP1, 0, 40, 90)
    translate(upperright_verticesP1, [+1+translationVectorPiglet[0],0+translationVectorPiglet[1],-0.85+translationVectorPiglet[2]])

    var vertex_bufferP8 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP8);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(upperright_verticesP1), GL.STATIC_DRAW);

    var color_bufferP8 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP8);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(upperright_colorP1), GL.STATIC_DRAW);

    var index_bufferP8 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP8);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(upperright_indicesP1), GL.STATIC_DRAW);

    // TANGAN KANAN BAWAH
    var icosahedronP9 = new Icosahedron3D(1);
    bottomright_verticesP1 = icosahedronP9.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    bottomright_indicesP1 = icosahedronP9.TriangleIndices;
    bottomright_colorP1 = generateColors(bottomright_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(bottomright_verticesP1, [0.5,0.35,1.1])
    rotate(bottomright_verticesP1, 50, -40, 90)
    translate(bottomright_verticesP1, [1+translationVectorPiglet[0],0.5+translationVectorPiglet[1],0.4+translationVectorPiglet[2]])
    
    var vertex_bufferP9 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP9);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bottomright_verticesP1), GL.STATIC_DRAW);
    
    var color_bufferP9 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP9);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bottomright_colorP1), GL.STATIC_DRAW);
    
    var index_bufferP9 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP9);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(bottomright_indicesP1), GL.STATIC_DRAW);
    
    // KAKI KIRI
    var L_foot = generateCube(-0.15, 0, 1, [0.45,0.3,0.7], normalizeColor(arrColors2[14]), false);
    scale(L_foot.vertices, [-2,2,2])
    rotate(L_foot.vertices,-10,0,0)
    translate(L_foot.vertices, [0+translationVectorPiglet[0],0+translationVectorPiglet[1],0+translationVectorPiglet[2]])
    var L_FOOT_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, L_FOOT_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(L_foot.vertices), GL.STATIC_DRAW);
    
    var L_FOOT_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, L_FOOT_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(L_foot.colors), GL.STATIC_DRAW);

    var L_FOOT_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, L_FOOT_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(L_foot.faces), GL.STATIC_DRAW);

    // KAKI KIRI BAWAH
    var icosahedronP10 = new Icosahedron3D(1);
    f_bottom_left_verticesP1 = icosahedronP10.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    f_bottom_left_indicesP1 = icosahedronP10.TriangleIndices;
    f_bottom_left_colorP1 = generateColors(f_bottom_left_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(f_bottom_left_verticesP1, [0.3,0.35,0.575])
    rotate(f_bottom_left_verticesP1, 0, 70, 90)
    translate(f_bottom_left_verticesP1, [-0.7+translationVectorPiglet[0],0+translationVectorPiglet[1],2.6+translationVectorPiglet[2]])

    var vertex_bufferP10 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP10);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(f_bottom_left_verticesP1), GL.STATIC_DRAW);

    var color_bufferP10 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP10);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(f_bottom_left_colorP1), GL.STATIC_DRAW);

    var index_bufferP10 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP10);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(f_bottom_left_indicesP1), GL.STATIC_DRAW);

    // KAKI KANAN
    var R_foot = generateCube(0.15, 0, 1, [0.45,0.3,0.7], normalizeColor(arrColors2[14]), false);
    scale(R_foot.vertices, [-2,2,2])
    rotate(R_foot.vertices,10,0,0)
    translate(R_foot.vertices, [0+translationVectorPiglet[0],0+translationVectorPiglet[1],0+translationVectorPiglet[2]])
    var R_FOOT_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, R_FOOT_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(R_foot.vertices), GL.STATIC_DRAW);
    
    var R_FOOT_COLORS = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, R_FOOT_COLORS);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(R_foot.colors), GL.STATIC_DRAW);

    var R_FOOT_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, R_FOOT_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(R_foot.faces), GL.STATIC_DRAW);

    // KAKI KANAN BAWAH
    var icosahedronP11 = new Icosahedron3D(1);
    f_bottom_right_verticesP1 = icosahedronP11.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    f_bottom_right_indicesP1 = icosahedronP11.TriangleIndices;
    f_bottom_right_colorP1 = generateColors(f_bottom_right_verticesP1, normalizeColor(arrColors2[14]), false);

    scale(f_bottom_right_verticesP1, [0.3,0.35,0.575])
    rotate(f_bottom_right_verticesP1, 0, -70, 90)
    translate(f_bottom_right_verticesP1, [0.7+translationVectorPiglet[0],0+translationVectorPiglet[1],2.6+translationVectorPiglet[2]])

    var vertex_bufferP11 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP11);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(f_bottom_right_verticesP1), GL.STATIC_DRAW);

    var color_bufferP11 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP11);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(f_bottom_right_colorP1), GL.STATIC_DRAW);

    var index_bufferP11 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP11);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(f_bottom_right_indicesP1), GL.STATIC_DRAW);

    // ALIS KANAN
    
        var alisKanan = generateQuadraticBezierCurve(controlPoints, .01);
        rotate(alisKanan.vertices, -30, 80, 45)
        scale(alisKanan.vertices, [.05,.1,.005])
        translate(alisKanan.vertices, [5+translationVectorPiglet[0],1.2+translationVectorPiglet[1], -2.7+translationVectorPiglet[2]])
        var vertex_alisKanan_buffer = GL.createBuffer();
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alisKanan_buffer);
        GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alisKanan.vertices), GL.STATIC_DRAW);
        var color_alisKanan_buffer = GL.createBuffer();
        GL.bindBuffer(GL.ARRAY_BUFFER, color_alisKanan_buffer);
        GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(generateColors(alisKanan.vertices, [0,0,0], false)), GL.STATIC_DRAW);

         // ALIS KIRI
        const controlPoints2 = [
        // Top left corner (A)
        [5, 12],
      
        // Top control points (additional)
        [6, 11.5],
        [7, 11],
      
        // Top middle control point (A1) - Adjust y-coordinate for more curvature
        [8, 10],
      
        // Top control points (additional)
        [9, 11],
        [10, 11.5],
      
        // Top right corner (B)
        [15, 12],
      
        // Bottom right corner (C)
        [15, 7],
      
        // Bottom control points (additional)
        [14, 7.5],
        [13, 8],
      
        // Bottom middle control point (C1) - Adjust y-coordinate for more curvature
        [8, 6],
      
        // Bottom control points (additional)
        [7, 8],
        [6, 7.5],
      
        // Bottom left corner (D)
        [5, 7],
      ];
      
    
        var alisKiri = generateQuadraticBezierCurve(controlPoints2, .01);
        rotate(alisKiri.vertices, -30, 80, 45)
        scale(alisKiri.vertices, [.04,.1,.005])
        translate(alisKiri.vertices, [4.3+translationVectorPiglet[0],1.2+translationVectorPiglet[1], -2.7+translationVectorPiglet[2]])
        var vertex_alisKiri_buffer = GL.createBuffer();
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alisKiri_buffer);
        GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alisKiri.vertices), GL.STATIC_DRAW);
        var color_alisKiri_buffer = GL.createBuffer();
        GL.bindBuffer(GL.ARRAY_BUFFER, color_alisKiri_buffer);
        GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(generateColors(alisKiri.vertices, [0,0,0], false)), GL.STATIC_DRAW);

    // HEAD
    var icosahedron = new Icosahedron3D(1);
    head_vertices = icosahedron.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    head_indices = icosahedron.TriangleIndices;
    head_color = generateColors(head_vertices, normalizeColor(arrColors[9]), false);

    // scale(head_vertices, [0,0,2])
    translate(head_vertices, [0,0,0], vecTransPooh)
    // rotate(head_vertices, 0, 0, 90)
    var vertex_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_vertices), GL.STATIC_DRAW);
    var color_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(head_color), GL.STATIC_DRAW);
    var index_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(head_indices), GL.STATIC_DRAW);
    
    // Jaw
    var icosahedron_jaw = new Icosahedron3D(1);
    jaw_vertices = icosahedron_jaw.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    jaw_indices = icosahedron_jaw.TriangleIndices;
    jaw_color = generateColors(jaw_vertices, normalizeColor(arrColors[9]), false);

    scale(jaw_vertices, [1.2,0,0])
    translate(jaw_vertices, [0,0,0.5], vecTransPooh)
    // rotate(jaw_vertices, 0, 0, 90)
    var vertex_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_jaw_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(jaw_vertices), GL.STATIC_DRAW);
    var color_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_jaw_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(jaw_color), GL.STATIC_DRAW);
    var index_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_jaw_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(jaw_indices), GL.STATIC_DRAW);
    
    // Congor
    var icosahedron_congor = new Icosahedron3D(1);
    congor_vertices = icosahedron_congor.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    congor_indices = icosahedron_congor.TriangleIndices;
    congor_color = generateColors(congor_vertices, normalizeColor(arrColors[10]), false);

    scale(congor_vertices, [0.7,.9,0.7])
    translate(congor_vertices, [0,0.5,0.5], vecTransPooh)
    // rotate(congor_vertices, 0, 0, 90)
    var vertex_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_congor_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(congor_vertices), GL.STATIC_DRAW);
    var color_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_congor_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(congor_color), GL.STATIC_DRAW);
    var index_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_congor_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(congor_indices), GL.STATIC_DRAW);
    
    // HIDUNG
    var icosahedron_hidung = new Icosahedron3D(1);
    hidung_vertices = icosahedron_hidung.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    hidung_indices = icosahedron_hidung.TriangleIndices;
    hidung_color = generateColors(hidung_vertices, normalizeColor(arrColors[8]));

    scale(hidung_vertices, [0.3,.2,0.2])
    translate(hidung_vertices, [0,1.2,.2], vecTransPooh)
    // rotate(hidung_vertices, 0, 0, 90)
    var vertex_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_hidung_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(hidung_vertices), GL.STATIC_DRAW);
    var color_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_hidung_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(hidung_color), GL.STATIC_DRAW);
    var index_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_hidung_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(hidung_indices), GL.STATIC_DRAW);
    
    // Telinga KIRI
    var icosahedron_telinga1 = new Icosahedron3D(1);
    telinga1_vertices = icosahedron_telinga1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    telinga1_indices = icosahedron_telinga1.TriangleIndices;
    telinga1_color = generateColors(telinga1_vertices, normalizeColor(arrColors[10]));

    scale(telinga1_vertices, [.3,.1,.2])
    rotate(telinga1_vertices, -60, 0, 0)
    translate(telinga1_vertices, [-.7,0,-1], vecTransPooh)
    var vertex_telinga1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_telinga1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(telinga1_vertices), GL.STATIC_DRAW);
    var color_telinga1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_telinga1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(telinga1_color), GL.STATIC_DRAW);
    var index_telinga1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_telinga1_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(telinga1_indices), GL.STATIC_DRAW);
    
    // Telinga KANAN
    var icosahedron_telinga2 = new Icosahedron3D(1);
    telinga2_vertices = icosahedron_telinga2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    telinga2_indices = icosahedron_telinga2.TriangleIndices;
    telinga2_color = generateColors(telinga2_vertices, normalizeColor(arrColors[10]));

    scale(telinga2_vertices, [.3,.1,.2])
    rotate(telinga2_vertices, 45, 0, 0)
    translate(telinga2_vertices, [.7,0,-1], vecTransPooh)
    var vertex_telinga2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_telinga2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(telinga2_vertices), GL.STATIC_DRAW);
    var color_telinga2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_telinga2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(telinga2_color), GL.STATIC_DRAW);
    var index_telinga2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_telinga2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(telinga2_indices), GL.STATIC_DRAW);
    
    // Mata kiri
    var icosahedron_mata1 = new Icosahedron3D(1);
    mata1_vertices = icosahedron_mata1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    mata1_indices = icosahedron_mata1.TriangleIndices;
    mata1_color = generateColors(mata1_vertices, arrColors[8]);

    scale(mata1_vertices, [.1,.1,.15])
    translate(mata1_vertices, [-.4,.7, -.6], vecTransPooh)
    // rotate(mata1_vertices, 0, 0, 90)
    var vertex_mata1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mata1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mata1_vertices), GL.STATIC_DRAW);
    var color_mata1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_mata1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mata1_color), GL.STATIC_DRAW);
    var index_mata1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_mata1_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(mata1_indices), GL.STATIC_DRAW);
    
    // Mata Kanan
    var icosahedron_mata2 = new Icosahedron3D(1);
    mata2_vertices = icosahedron_mata2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    mata2_indices = icosahedron_mata2.TriangleIndices;
    mata2_color = generateColors(mata2_vertices, arrColors[8]);

    scale(mata2_vertices, [.1,.1,.15])
    translate(mata2_vertices, [.4,.7, -.6], vecTransPooh)
    // rotate(mata2_vertices, 0, 0, 90)
    var vertex_mata2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mata2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mata2_vertices), GL.STATIC_DRAW);
    var color_mata2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_mata2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mata2_color), GL.STATIC_DRAW);
    var index_mata2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_mata2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(mata2_indices), GL.STATIC_DRAW);
    
    // Alis Kanan
    var icosahedron_alis1 = new Icosahedron3D(1);
    alis1_vertices = icosahedron_alis1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    alis1_indices = icosahedron_alis1.TriangleIndices;
    alis1_color = generateColors(alis1_vertices, arrColors[8]);

    scale(alis1_vertices, [.2,.01,.02])
    rotate(alis1_vertices, -20, 0, 0)
    translate(alis1_vertices, [.45,.5, -.75], vecTransPooh)
    var vertex_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alis1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alis1_vertices), GL.STATIC_DRAW);
    var color_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_alis1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alis1_color), GL.STATIC_DRAW);
    var index_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_alis1_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(alis1_indices), GL.STATIC_DRAW);
    
    // Alis Kiri
    var icosahedron_alis2 = new Icosahedron3D(1);
    alis2_vertices = icosahedron_alis2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    alis2_indices = icosahedron_alis2.TriangleIndices;
    alis2_color = generateColors(alis2_vertices, arrColors[8]);

    scale(alis2_vertices, [.2,.01,.02])
    rotate(alis2_vertices, 30, 0, 0)
    translate(alis2_vertices, [-.45,.5, -.75], vecTransPooh)
    var vertex_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alis2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alis2_vertices), GL.STATIC_DRAW);
    var color_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_alis2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(alis2_color), GL.STATIC_DRAW);
    var index_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_alis2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(alis2_indices), GL.STATIC_DRAW);
    
    // Badan Pooh
    var icosahedron_badan = new Icosahedron3D(1);
    badan_vertices = icosahedron_badan.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    badan_indices = icosahedron_badan.TriangleIndices;
    badan_color = generateColors(badan_vertices, normalizeColor(arrColors[9]));

    scale(badan_vertices, [1.5,1,2])
    // rotate(badan_vertices, 30, 0, 0)
    translate(badan_vertices, [0,0, 2.5], vecTransPooh)
    var vertex_badan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_vertices), GL.STATIC_DRAW);
    var color_badan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_color), GL.STATIC_DRAW);
    var index_badan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(badan_indices), GL.STATIC_DRAW);
    
    // Tangan Atas Kiri
    var poohtak = new Cylinder3D(.3,2,100,100, true);
    poohtak_vertices = poohtak.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohtak_indices = poohtak.TriangleIndices;
    poohtak_color = generateColors(poohtak_vertices, normalizeColor(arrColors[9]));
    // scale(poohtak_vertices, [1.5,1,2])
    rotate(poohtak_vertices, -45, 90, 0)
    translate(poohtak_vertices, [-1.5,0, 2], vecTransPooh)
    var vertex_poohtak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtak_vertices), GL.STATIC_DRAW);
    
    var color_poohtak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtak_color), GL.STATIC_DRAW);

    var index_poohtak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtak_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohtak_indices), GL.STATIC_DRAW);
    
    // Tangan Atas Kanan
    var poohtakn = new Cylinder3D(.3,2,100,100, true);
    poohtakn_vertices = poohtakn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohtakn_indices = poohtakn.TriangleIndices;
    poohtakn_color = generateColors(poohtakn_vertices, normalizeColor(arrColors[9]));
    // scale(poohtakn_vertices, [1.5,1,2])
    rotate(poohtakn_vertices, 120, 90, 0)
    translate(poohtakn_vertices, [1.5,0, 1], vecTransPooh)
    var vertex_poohtakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtakn_vertices), GL.STATIC_DRAW);
    
    var color_poohtakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtakn_color), GL.STATIC_DRAW);

    var index_poohtakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtakn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohtakn_indices), GL.STATIC_DRAW);
    
    // Kaki Atas Kiri
    var poohkak = new Cylinder3D(.3,2,100,100, true);
    poohkak_vertices = poohkak.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohkak_indices = poohkak.TriangleIndices;
    poohkak_color = generateColors(poohkak_vertices, normalizeColor(arrColors[9]));
    // scale(poohkak_vertices, [1.5,1,2])
    rotate(poohkak_vertices, -60, 90, -90)
    translate(poohkak_vertices, [-.8,.5, 4], vecTransPooh)
    var vertex_poohkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkak_vertices), GL.STATIC_DRAW);
    
    var color_poohkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkak_color), GL.STATIC_DRAW);

    var index_poohkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkak_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohkak_indices), GL.STATIC_DRAW);
    
    // Kaki Atas Kanan
    var poohkakn = new Cylinder3D(.3,2,100,100, true);
    poohkakn_vertices = poohkakn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohkakn_indices = poohkakn.TriangleIndices;
    poohkakn_color = generateColors(poohkakn_vertices, normalizeColor(arrColors[9]));
    // scale(poohkakn_vertices, [1.5,1,2])
    rotate(poohkakn_vertices, -60, 90, -90)
    translate(poohkakn_vertices, [.8,.5, 4], vecTransPooh)
    var vertex_poohkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkakn_vertices), GL.STATIC_DRAW);
    
    var color_poohkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkakn_color), GL.STATIC_DRAW);

    var index_poohkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkakn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohkakn_indices), GL.STATIC_DRAW);
    
    // Tangan Bawah Kiri
    var poohtbk = new Icosahedron3D(1);
    poohtbk_vertices = poohtbk.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohtbk_indices = poohtbk.TriangleIndices;
    poohtbk_color = generateColors(poohtbk_vertices, normalizeColor(arrColors[10]));
    scale(poohtbk_vertices, [.35,.35,.35])
    // rotate(poohtbk_vertices, -60, 90, -90)
    translate(poohtbk_vertices, [2.5,0, .4], vecTransPooh)
    var vertex_poohtbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtbk_vertices), GL.STATIC_DRAW);
    
    var color_poohtbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtbk_color), GL.STATIC_DRAW);

    var index_poohtbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtbk_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohtbk_indices), GL.STATIC_DRAW);
    
    // Tangan Bawah Kanan
    var poohtbkn = new Icosahedron3D(1);
    poohtbkn_vertices = poohtbkn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohtbkn_indices = poohtbkn.TriangleIndices;
    poohtbkn_color = generateColors(poohtbkn_vertices, normalizeColor(arrColors[10]));
    scale(poohtbkn_vertices, [.35,.35,.35])
    // rotate(poohtbkn_vertices, -60, 90, -90)
    translate(poohtbkn_vertices, [-2.3,0, 2.8], vecTransPooh)
    var vertex_poohtbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtbkn_vertices), GL.STATIC_DRAW);
    
    var color_poohtbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohtbkn_color), GL.STATIC_DRAW);

    var index_poohtbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtbkn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohtbkn_indices), GL.STATIC_DRAW);
    
    // kbk
    var poohkbk = new Icosahedron3D(1);
    poohkbk_vertices = poohkbk.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohkbk_indices = poohkbk.TriangleIndices;
    poohkbk_color = generateColors(poohkbk_vertices, normalizeColor(arrColors[10]));
    scale(poohkbk_vertices, [.35,.5,.35])
    // rotate(poohkbk_vertices, -60, 90, -90)
    translate(poohkbk_vertices, [-.8,1.3, 4.5], vecTransPooh)
    var vertex_poohkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkbk_vertices), GL.STATIC_DRAW);
    
    var color_poohkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkbk_color), GL.STATIC_DRAW);

    var index_poohkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkbk_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohkbk_indices), GL.STATIC_DRAW);
    
    // kbkn
    var poohkbkn = new Icosahedron3D(1);
    poohkbkn_vertices = poohkbkn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    poohkbkn_indices = poohkbkn.TriangleIndices;
    poohkbkn_color = generateColors(poohkbkn_vertices, normalizeColor(arrColors[10]));
    scale(poohkbkn_vertices, [.35,.5,.35])
    // rotate(poohkbkn_vertices, -60, 90, -90)
    translate(poohkbkn_vertices, [.8,1.3, 4.5], vecTransPooh)
    var vertex_poohkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkbkn_vertices), GL.STATIC_DRAW);
    
    var color_poohkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(poohkbkn_color), GL.STATIC_DRAW);

    var index_poohkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkbkn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(poohkbkn_indices), GL.STATIC_DRAW);

    // Mulut
    // Control points for smiling Pooh mouth (adjust as needed)
    
  

    var mulut_pooh = generateQuadraticBezierCurve(controlPoints, .01);
    rotate(mulut_pooh.vertices, -30, 80, 45)
    scale(mulut_pooh.vertices, [.07,.07,.07])
    translate(mulut_pooh.vertices, [-.4,1.63, -.3], vecTransPooh)
    var vertex_mulut_pooh_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_pooh_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mulut_pooh.vertices), GL.STATIC_DRAW);
    var color_mulut_pooh_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_pooh_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(generateColors(mulut_pooh.vertices, [0,0,0], false)), GL.STATIC_DRAW);
    
    // Mulut
    // Control points for smiling Pooh mouth (adjust as needed)
    
  

    var mulut_piglet = generateQuadraticBezierCurve(controlPoints, .01);
    rotate(mulut_piglet.vertices, -30, 80, 45)
    scale(mulut_piglet.vertices, [.07,.07,.07])
    translate(mulut_piglet.vertices, [-.4,1, -3], translationVectorPiglet)
    var vertex_mulut_piglet_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_piglet_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mulut_piglet.vertices), GL.STATIC_DRAW);
    var color_mulut_piglet_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_piglet_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(generateColors(mulut_piglet.vertices, [0,0,0], false)), GL.STATIC_DRAW);


    // RABBIT
    // HEAD
    var rabbit = new Icosahedron3D(1);
    rabbit_head_vertices = rabbit.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_head_indices = rabbit.TriangleIndices;
    rabbit_head_color = generateColors(rabbit_head_vertices, normalizeColor(arrColors[14]), false);

    // scale(rabbit_head_vertices, [0,0,2])
    translate(rabbit_head_vertices, [0,0,0], vecTransRabbit)
    // rotate(rabbit_head_vertices, 0, 0, 90)
    var vertex_rabbit_head_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_head_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_head_vertices), GL.STATIC_DRAW);
    var color_rabbit_head_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_head_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_head_color), GL.STATIC_DRAW);
    var index_rabbit_head_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_head_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_head_indices), GL.STATIC_DRAW);
    
    // Jaw
    var rabbit_jaw = new Icosahedron3D(1);
    rabbit_jaw_vertices = rabbit_jaw.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_jaw_indices = rabbit_jaw.TriangleIndices;
    rabbit_jaw_color = generateColors(rabbit_jaw_vertices, normalizeColor(arrColors[14]), false);

    scale(rabbit_jaw_vertices, [1.2,0,0])
    translate(rabbit_jaw_vertices, [0,0,0.5], vecTransRabbit)
    // rotate(rabbit_jaw_vertices, 0, 0, 90)
    var vertex_rabbit_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_jaw_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_jaw_vertices), GL.STATIC_DRAW);
    var color_rabbit_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_jaw_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_jaw_color), GL.STATIC_DRAW);
    var index_rabbit_jaw_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_jaw_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_jaw_indices), GL.STATIC_DRAW);
    
    // Congor
    var rabbit_congor = new Icosahedron3D(1);
    rabbit_congor_vertices = rabbit_congor.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_congor_indices = rabbit_congor.TriangleIndices;
    rabbit_congor_color = generateColors(rabbit_congor_vertices, (arrColors[6]), false);

    scale(rabbit_congor_vertices, [0.7,.9,0.7])
    translate(rabbit_congor_vertices, [0,0.5,0.5], vecTransRabbit)
    // rotate(rabbit_congor_vertices, 0, 0, 90)
    var vertex_rabbit_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_congor_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_congor_vertices), GL.STATIC_DRAW);
    var color_rabbit_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_congor_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_congor_color), GL.STATIC_DRAW);
    var index_rabbit_congor_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_congor_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_congor_indices), GL.STATIC_DRAW);
    
    // HIDUNG
    var rabbit_hidung = new Icosahedron3D(1);
    rabbit_hidung_vertices = rabbit_hidung.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_hidung_indices = rabbit_hidung.TriangleIndices;
    rabbit_hidung_color = generateColors(rabbit_hidung_vertices, normalizeColor(arrColors[15]));

    scale(rabbit_hidung_vertices, [0.3,.2,0.2])
    translate(rabbit_hidung_vertices, [0,1.2,.2], vecTransRabbit)
    // rotate(rabbit_hidung_vertices, 0, 0, 90)
    var vertex_rabbit_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_hidung_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_hidung_vertices), GL.STATIC_DRAW);
    var color_rabbit_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_hidung_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_hidung_color), GL.STATIC_DRAW);
    var index_rabbit_hidung_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_hidung_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_hidung_indices), GL.STATIC_DRAW);
    
    // ML_rabbit
    var rabbit_mata_luar = new Icosahedron3D(1);
    rabbit_mata_luar_vertices = rabbit_mata_luar.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_mata_luar_indices = rabbit_mata_luar.TriangleIndices;
    rabbit_mata_luar_color = generateColors(rabbit_mata_luar_vertices, (arrColors[6]));

    scale(rabbit_mata_luar_vertices, [0.3,.3,0.35])
    translate(rabbit_mata_luar_vertices, [-.3,.6,-.3], vecTransRabbit)
    // rotate(rabbit_mata_luar_vertices, 0, 0, 90)
    var vertex_rabbit_mata_luar_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_luar_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_luar_vertices), GL.STATIC_DRAW);
    var color_rabbit_mata_luar_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_luar_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_luar_color), GL.STATIC_DRAW);
    var index_rabbit_mata_luar_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_luar_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_mata_luar_indices), GL.STATIC_DRAW);
   
    // MD_rabbit
    var rabbit_mata_dalam = new Icosahedron3D(1);
    rabbit_mata_dalam_vertices = rabbit_mata_dalam.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_mata_dalam_indices = rabbit_mata_dalam.TriangleIndices;
    rabbit_mata_dalam_color = generateColors(rabbit_mata_dalam_vertices, normalizeColor(arrColors[8]));

    scale(rabbit_mata_dalam_vertices, [0.1,.1,0.1])
    translate(rabbit_mata_dalam_vertices, [-.35,.75,-.5], vecTransRabbit)
    // rotate(rabbit_mata_dalam_vertices, 0, 0, 90)
    var vertex_rabbit_mata_dalam_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_dalam_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_dalam_vertices), GL.STATIC_DRAW);
    var color_rabbit_mata_dalam_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_dalam_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_dalam_color), GL.STATIC_DRAW);
    var index_rabbit_mata_dalam_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_dalam_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_mata_dalam_indices), GL.STATIC_DRAW);
    
    // ML_rabbit KN
    var rabbit_mata_luar_kanan = new Icosahedron3D(1);
    rabbit_mata_luar_kanan_vertices = rabbit_mata_luar_kanan.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_mata_luar_kanan_indices = rabbit_mata_luar_kanan.TriangleIndices;
    rabbit_mata_luar_kanan_color = generateColors(rabbit_mata_luar_kanan_vertices, (arrColors[6]));

    scale(rabbit_mata_luar_kanan_vertices, [0.3,.3,0.35])
    translate(rabbit_mata_luar_kanan_vertices, [.3,.6,-.3], vecTransRabbit)
    // rotate(rabbit_mata_luar_kanan_vertices, 0, 0, 90)
    var vertex_rabbit_mata_luar_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_luar_kanan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_luar_kanan_vertices), GL.STATIC_DRAW);
    var color_rabbit_mata_luar_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_luar_kanan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_luar_kanan_color), GL.STATIC_DRAW);
    var index_rabbit_mata_luar_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_luar_kanan_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_mata_luar_kanan_indices), GL.STATIC_DRAW);
   
    // MD_rabbit KN
    var rabbit_mata_dalam_kanan = new Icosahedron3D(1);
    rabbit_mata_dalam_kanan_vertices = rabbit_mata_dalam_kanan.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_mata_dalam_kanan_indices = rabbit_mata_dalam_kanan.TriangleIndices;
    rabbit_mata_dalam_kanan_color = generateColors(rabbit_mata_dalam_kanan_vertices, normalizeColor(arrColors[8]));

    scale(rabbit_mata_dalam_kanan_vertices, [0.1,.1,0.1])
    translate(rabbit_mata_dalam_kanan_vertices, [.35,.75,-.5], vecTransRabbit)
    // rotate(rabbit_mata_dalam_kanan_vertices, 0, 0, 90)
    var vertex_rabbit_mata_dalam_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_dalam_kanan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_dalam_kanan_vertices), GL.STATIC_DRAW);
    var color_rabbit_mata_dalam_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_dalam_kanan_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_mata_dalam_kanan_color), GL.STATIC_DRAW);
    var index_rabbit_mata_dalam_kanan_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_dalam_kanan_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_mata_dalam_kanan_indices), GL.STATIC_DRAW);
    
    // rabbit_telinga
    var rabbit_telinga = new Icosahedron3D(1);
    rabbit_telinga_vertices = rabbit_telinga.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_telinga_indices = rabbit_telinga.TriangleIndices;
    rabbit_telinga_color = generateColors(rabbit_telinga_vertices, normalizeColor(arrColors[14]));

    scale(rabbit_telinga_vertices, [0.3,.1,1])
    rotate(rabbit_telinga_vertices, 30, 0, 0)
    translate(rabbit_telinga_vertices, [-.7,.2,-1.5], vecTransRabbit)
    var vertex_rabbit_telinga_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_vertices), GL.STATIC_DRAW);
    var color_rabbit_telinga_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_color), GL.STATIC_DRAW);
    var index_rabbit_telinga_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_telinga_indices), GL.STATIC_DRAW);
    
    // rabbit_telinga_dl
    var rabbit_telinga_dl = new Icosahedron3D(1);
    rabbit_telinga_dl_vertices = rabbit_telinga_dl.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_telinga_dl_indices = rabbit_telinga_dl.TriangleIndices;
    rabbit_telinga_dl_color = generateColors(rabbit_telinga_dl_vertices, normalizeColor(arrColors[15]));

    scale(rabbit_telinga_dl_vertices, [0.2,.05,.8])
    rotate(rabbit_telinga_dl_vertices, 30, 0, 0)
    translate(rabbit_telinga_dl_vertices, [-.7,.28,-1.5], vecTransRabbit)
    var vertex_rabbit_telinga_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_dl_vertices), GL.STATIC_DRAW);
    var color_rabbit_telinga_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_dl_color), GL.STATIC_DRAW);
    var index_rabbit_telinga_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_dl_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_telinga_dl_indices), GL.STATIC_DRAW);
    
    // rabbit_telinga_kn
    var rabbit_telinga_kn = new Icosahedron3D(1);
    rabbit_telinga_kn_vertices = rabbit_telinga_kn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_telinga_kn_indices = rabbit_telinga_kn.TriangleIndices;
    rabbit_telinga_kn_color = generateColors(rabbit_telinga_kn_vertices, normalizeColor(arrColors[14]));

    scale(rabbit_telinga_kn_vertices, [0.3,.1,1])
    rotate(rabbit_telinga_kn_vertices, -30, 0, 0)
    translate(rabbit_telinga_kn_vertices, [.7,.2,-1.5], vecTransRabbit)
    var vertex_rabbit_telinga_kn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_kn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_kn_vertices), GL.STATIC_DRAW);
    var color_rabbit_telinga_kn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_kn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_kn_color), GL.STATIC_DRAW);
    var index_rabbit_telinga_kn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_kn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_telinga_kn_indices), GL.STATIC_DRAW);
    
    // rabbit_telinga_kn_dl
    var rabbit_telinga_kn_dl = new Icosahedron3D(1);
    rabbit_telinga_kn_dl_vertices = rabbit_telinga_kn_dl.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_telinga_kn_dl_indices = rabbit_telinga_kn_dl.TriangleIndices;
    rabbit_telinga_kn_dl_color = generateColors(rabbit_telinga_kn_dl_vertices, normalizeColor(arrColors[15]));

    scale(rabbit_telinga_kn_dl_vertices, [0.2,.05,.8])
    rotate(rabbit_telinga_kn_dl_vertices, -30, 0, 0)
    translate(rabbit_telinga_kn_dl_vertices, [.7,.28,-1.5], vecTransRabbit)
    var vertex_rabbit_telinga_kn_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_kn_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_kn_dl_vertices), GL.STATIC_DRAW);
    var color_rabbit_telinga_kn_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_kn_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_telinga_kn_dl_color), GL.STATIC_DRAW);
    var index_rabbit_telinga_kn_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_kn_dl_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_telinga_kn_dl_indices), GL.STATIC_DRAW);

    // Alis Kanan
    var icosahedron_rabbit_alis1 = new Icosahedron3D(1);
    rabbit_alis1_vertices = icosahedron_rabbit_alis1.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_alis1_indices = icosahedron_rabbit_alis1.TriangleIndices;
    rabbit_alis1_color = generateColors(rabbit_alis1_vertices, arrColors[8]);

    scale(rabbit_alis1_vertices, [.2,.01,.02])
    rotate(rabbit_alis1_vertices, -20, 0, 0)
    translate(rabbit_alis1_vertices, [.45,.5, -.75], vecTransRabbit)
    var vertex_rabbit_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_alis1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_alis1_vertices), GL.STATIC_DRAW);
    var color_rabbit_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_alis1_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_alis1_color), GL.STATIC_DRAW);
    var index_rabbit_alis1_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_alis1_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_alis1_indices), GL.STATIC_DRAW);
    
    // Alis Kiri
    var icosahedron_rabbit_alis2 = new Icosahedron3D(1);
    rabbit_alis2_vertices = icosahedron_rabbit_alis2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbit_alis2_indices = icosahedron_rabbit_alis2.TriangleIndices;
    rabbit_alis2_color = generateColors(rabbit_alis2_vertices, arrColors[8]);

    scale(rabbit_alis2_vertices, [.2,.01,.02])
    rotate(rabbit_alis2_vertices, 30, 0, 0)
    translate(rabbit_alis2_vertices, [-.45,.5, -.75], vecTransRabbit)
    var vertex_rabbit_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_alis2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_alis2_vertices), GL.STATIC_DRAW);
    var color_rabbit_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_alis2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbit_alis2_color), GL.STATIC_DRAW);
    var index_rabbit_alis2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_alis2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbit_alis2_indices), GL.STATIC_DRAW);

    // mulut_rabbit
    var mulut_rabbit = generateQuadraticBezierCurve(controlPoints, .01);
    rotate(mulut_rabbit.vertices, -30, 80, 45)
    scale(mulut_rabbit.vertices, [.07,.07,.07])
    translate(mulut_rabbit.vertices, [-.4,1.63, -.3], vecTransRabbit)
    var vertex_mulut_rabbit_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_rabbit_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mulut_rabbit.vertices), GL.STATIC_DRAW);
    var color_mulut_rabbit_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_rabbit_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(generateColors(mulut_rabbit.vertices, [0,0,0], false)), GL.STATIC_DRAW);

    // badan_rabbit
    var icosahedron_badan_rabbit = new Icosahedron3D(1);
    badan_rabbit_vertices = icosahedron_badan_rabbit.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    badan_rabbit_indices = icosahedron_badan_rabbit.TriangleIndices;
    badan_rabbit_color = generateColors(badan_rabbit_vertices, normalizeColor(arrColors[14]));

    scale(badan_rabbit_vertices, [1,1,2])
    // rotate(badan_rabbit_vertices, 30, 0, 0)
    translate(badan_rabbit_vertices, [0,0, 2.5], vecTransRabbit)
    var vertex_badan_rabbit_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_vertices), GL.STATIC_DRAW);
    var color_badan_rabbit_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_color), GL.STATIC_DRAW);
    var index_badan_rabbit_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(badan_rabbit_indices), GL.STATIC_DRAW);
    
    // badan_rabbit_2
    var icosahedron_badan_rabbit2 = new Icosahedron3D(1);
    badan_rabbit2_vertices = icosahedron_badan_rabbit2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    badan_rabbit2_indices = icosahedron_badan_rabbit2.TriangleIndices;
    badan_rabbit2_color = generateColors(badan_rabbit2_vertices, normalizeColor(arrColors[14]));

    scale(badan_rabbit2_vertices, [1.3,1,1.5])
    // rotate(badan_rabbit2_vertices, 30, 0, 0)
    translate(badan_rabbit2_vertices, [0,0, 3.5], vecTransRabbit)
    var vertex_badan_rabbit2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit2_vertices), GL.STATIC_DRAW);
    var color_badan_rabbit2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit2_color), GL.STATIC_DRAW);
    var index_badan_rabbit2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(badan_rabbit2_indices), GL.STATIC_DRAW);
    
    // badan_rabbit_dl
    var icosahedron_badan_rabbit_dl = new Icosahedron3D(1);
    badan_rabbit_dl_vertices = icosahedron_badan_rabbit_dl.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    badan_rabbit_dl_indices = icosahedron_badan_rabbit_dl.TriangleIndices;
    badan_rabbit_dl_color = generateColors(badan_rabbit_dl_vertices, (arrColors[14]));

    scale(badan_rabbit_dl_vertices, [.8,1,1.8])
    // rotate(badan_rabbit_dl_vertices, 30, 0, 0)
    translate(badan_rabbit_dl_vertices, [0,0.2, 2.5], vecTransRabbit)
    var vertex_badan_rabbit_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_dl_vertices), GL.STATIC_DRAW);
    var color_badan_rabbit_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_dl_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_dl_color), GL.STATIC_DRAW);
    var index_badan_rabbit_dl_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_dl_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(badan_rabbit_dl_indices), GL.STATIC_DRAW);
    
    // badan_rabbit_dl_2
    var icosahedron_badan_rabbit_dl2 = new Icosahedron3D(1);
    badan_rabbit_dl2_vertices = icosahedron_badan_rabbit_dl2.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    badan_rabbit_dl2_indices = icosahedron_badan_rabbit_dl2.TriangleIndices;
    badan_rabbit_dl2_color = generateColors(badan_rabbit_dl2_vertices, (arrColors[14]));

    scale(badan_rabbit_dl2_vertices, [1.1,1,1.5])
    // rotate(badan_rabbit_dl2_vertices, 30, 0, 0)
    translate(badan_rabbit_dl2_vertices, [0,0.2, 3.5], vecTransRabbit)
    var vertex_badan_rabbit_dl2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_dl2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_dl2_vertices), GL.STATIC_DRAW);
    var color_badan_rabbit_dl2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_dl2_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(badan_rabbit_dl2_color), GL.STATIC_DRAW);
    var index_badan_rabbit_dl2_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_dl2_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(badan_rabbit_dl2_indices), GL.STATIC_DRAW);
    
    // rabbittak
    var rabbittak = new Cylinder3D(.3,2,100,100, true);
    rabbittak_vertices = rabbittak.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbittak_indices = rabbittak.TriangleIndices;
    rabbittak_color = generateColors(rabbittak_vertices, normalizeColor(arrColors[14]));
    // scale(rabbittak_vertices, [1.5,1,2])
    rotate(rabbittak_vertices, 45, 90, 0)
    translate(rabbittak_vertices, [-1.5,0, 1], vecTransRabbit)
    var vertex_rabbittak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittak_vertices), GL.STATIC_DRAW);
    
    var color_rabbittak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittak_color), GL.STATIC_DRAW);

    var index_rabbittak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittak_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbittak_indices), GL.STATIC_DRAW);
    
    // Tangan Atas Kanan
    var rabbittakn = new Cylinder3D(.3,2,100,100, true);
    rabbittakn_vertices = rabbittakn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbittakn_indices = rabbittakn.TriangleIndices;
    rabbittakn_color = generateColors(rabbittakn_vertices, normalizeColor(arrColors[14]));
    // scale(rabbittakn_vertices, [1.5,1,2])
    rotate(rabbittakn_vertices, 120, 90, 0)
    translate(rabbittakn_vertices, [1.5,0, 1], vecTransRabbit)
    var vertex_rabbittakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittakn_vertices), GL.STATIC_DRAW);
    
    var color_rabbittakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittakn_color), GL.STATIC_DRAW);

    var index_rabbittakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittakn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbittakn_indices), GL.STATIC_DRAW);
    
    // Kaki Atas Kiri
    var rabbitkak = new Cylinder3D(.3,2,100,100, true);
    rabbitkak_vertices = rabbitkak.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbitkak_indices = rabbitkak.TriangleIndices;
    rabbitkak_color = generateColors(rabbitkak_vertices, normalizeColor(arrColors[14]));
    // scale(rabbitkak_vertices, [1.5,1,2])
    rotate(rabbitkak_vertices, -60, 90, -90)
    translate(rabbitkak_vertices, [-.8,.5, 4], vecTransRabbit)
    var vertex_rabbitkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkak_vertices), GL.STATIC_DRAW);
    
    var color_rabbitkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkak_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkak_color), GL.STATIC_DRAW);

    var index_rabbitkak_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkak_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbitkak_indices), GL.STATIC_DRAW);
    
    // Kaki Atas Kanan
    var rabbitkakn = new Cylinder3D(.3,2,100,100, true);
    rabbitkakn_vertices = rabbitkakn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbitkakn_indices = rabbitkakn.TriangleIndices;
    rabbitkakn_color = generateColors(rabbitkakn_vertices, normalizeColor(arrColors[14]));
    // scale(rabbitkakn_vertices, [1.5,1,2])
    rotate(rabbitkakn_vertices, -60, 90, -90)
    translate(rabbitkakn_vertices, [.8,.5, 4], vecTransRabbit)
    var vertex_rabbitkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkakn_vertices), GL.STATIC_DRAW);
    
    var color_rabbitkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkakn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkakn_color), GL.STATIC_DRAW);

    var index_rabbitkakn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkakn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbitkakn_indices), GL.STATIC_DRAW);
    
    // Tangan Bawah Kiri
    var rabbittbk = new Icosahedron3D(1);
    rabbittbk_vertices = rabbittbk.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbittbk_indices = rabbittbk.TriangleIndices;
    rabbittbk_color = generateColors(rabbittbk_vertices, normalizeColor(arrColors[15]));
    scale(rabbittbk_vertices, [.35,.35,.35])
    // rotate(rabbittbk_vertices, -60, 90, -90)
    translate(rabbittbk_vertices, [2.5,0, .4], vecTransRabbit)
    var vertex_rabbittbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittbk_vertices), GL.STATIC_DRAW);
    
    var color_rabbittbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittbk_color), GL.STATIC_DRAW);

    var index_rabbittbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittbk_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbittbk_indices), GL.STATIC_DRAW);
    
    // Tangan Bawah Kanan
    var rabbittbkn = new Icosahedron3D(1);
    rabbittbkn_vertices = rabbittbkn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbittbkn_indices = rabbittbkn.TriangleIndices;
    rabbittbkn_color = generateColors(rabbittbkn_vertices, normalizeColor(arrColors[15]));
    scale(rabbittbkn_vertices, [.35,.35,.35])
    // rotate(rabbittbkn_vertices, -60, 90, -90)
    translate(rabbittbkn_vertices, [-2.3,0, .25], vecTransRabbit)
    var vertex_rabbittbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittbkn_vertices), GL.STATIC_DRAW);
    
    var color_rabbittbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbittbkn_color), GL.STATIC_DRAW);

    var index_rabbittbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittbkn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbittbkn_indices), GL.STATIC_DRAW);
    
    // kbk
    var rabbitkbk = new Icosahedron3D(1);
    rabbitkbk_vertices = rabbitkbk.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbitkbk_indices = rabbitkbk.TriangleIndices;
    rabbitkbk_color = generateColors(rabbitkbk_vertices, normalizeColor(arrColors[15]));
    scale(rabbitkbk_vertices, [.35,.5,.35])
    // rotate(rabbitkbk_vertices, -60, 90, -90)
    translate(rabbitkbk_vertices, [-.8,1.3, 4.5], vecTransRabbit)
    var vertex_rabbitkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkbk_vertices), GL.STATIC_DRAW);
    
    var color_rabbitkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkbk_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkbk_color), GL.STATIC_DRAW);

    var index_rabbitkbk_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkbk_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbitkbk_indices), GL.STATIC_DRAW);
    
    // kbkn
    var rabbitkbkn = new Icosahedron3D(1);
    rabbitkbkn_vertices = rabbitkbkn.Points.reduce(function (a, b, i) { return i === 1 ? [a.x, a.y, a.z, b.x, b.y, b.z] : a.concat([b.x, b.y, b.z]); });
    rabbitkbkn_indices = rabbitkbkn.TriangleIndices;
    rabbitkbkn_color = generateColors(rabbitkbkn_vertices, normalizeColor(arrColors[15]));
    scale(rabbitkbkn_vertices, [.35,.5,.35])
    // rotate(rabbitkbkn_vertices, -60, 90, -90)
    translate(rabbitkbkn_vertices, [.8,1.3, 4.5], vecTransRabbit)
    var vertex_rabbitkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkbkn_vertices), GL.STATIC_DRAW);
    
    var color_rabbitkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkbkn_buffer);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rabbitkbkn_color), GL.STATIC_DRAW);

    var index_rabbitkbkn_buffer = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkbkn_buffer);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rabbitkbkn_indices), GL.STATIC_DRAW);


        // SUN
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferSUN);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferSUN);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferSUN);

        GL.drawElements(GL.TRIANGLES, sun_indices.length, GL.UNSIGNED_SHORT, 0);

        // CLOUD BULAT 1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB1);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB1);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB1);

        GL.drawElements(GL.TRIANGLES, cloudB1_indices.length, GL.UNSIGNED_SHORT, 0);

        // CLOUD BULAT 2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB2);

        GL.drawElements(GL.TRIANGLES, cloudB2_indices.length, GL.UNSIGNED_SHORT, 0);

        // CLOUD BULAT 3
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB3);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB3);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB3);

        GL.drawElements(GL.TRIANGLES, cloudB3_indices.length, GL.UNSIGNED_SHORT, 0);

        
        // CLOUD BULAT 4
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB4);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB4);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB4);

        GL.drawElements(GL.TRIANGLES, cloudB4_indices.length, GL.UNSIGNED_SHORT, 0);

        // CLOUD BULAT 5
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB5);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB5);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB5);

        GL.drawElements(GL.TRIANGLES, cloudB5_indices.length, GL.UNSIGNED_SHORT, 0);

        // CLOUD BULAT 6
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffercloudB6);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffercloudB6);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffercloudB6);

        GL.drawElements(GL.TRIANGLES, cloudB6_indices.length, GL.UNSIGNED_SHORT, 0);

        // SURFACE
        GL.bindBuffer(GL.ARRAY_BUFFER, SURFACE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, SURFACE_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, SURFACE_FACES);
        
        GL.drawElements(GL.TRIANGLES, objSurface.faces.length, GL.UNSIGNED_SHORT, 0);
        
        // TBRANCH
        GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TBRANCH_FACES);

        GL.drawElements(GL.TRIANGLES, objTb.faces.length, GL.UNSIGNED_SHORT, 0);

        // TBRANCH
        GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, TBRANCH2_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TBRANCH2_FACES);

        //GL.drawElements(GL.TRIANGLES, objTb2.faces.length, GL.UNSIGNED_SHORT, 0);
        
        // TLEAF
        GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TLEAF_FACES);

        GL.drawElements(GL.TRIANGLES, objLeaf.faces.length, GL.UNSIGNED_SHORT, 0);
        
        // TLEAF
        GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, TLEAF2_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TLEAF2_FACES);

        GL.drawElements(GL.TRIANGLES, objLeaf2.faces.length, GL.UNSIGNED_SHORT, 0);
        
        // HEAD PIGLET
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP1);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP1);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP1);

        GL.drawElements(GL.TRIANGLES, head_indicesP1.length, GL.UNSIGNED_SHORT, 0);

         // BODY PIGLET
         GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP2);
         GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
 
         GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP2);
         GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
 
         GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP2);
 
         GL.drawElements(GL.TRIANGLES, head_indicesP2.length, GL.UNSIGNED_SHORT, 0);
         
         // NOSE PIGLET
         GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX);
         GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
 
         GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLORS);
         GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
 
         GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_FACES);
 
         GL.drawElements(GL.TRIANGLES, nose.faces.length, GL.UNSIGNED_SHORT, 0);

         // NOSE PIGLET BULAT
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP5);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP5);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP5);

        GL.drawElements(GL.TRIANGLES, lefteye_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // MATA KIRI
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP3);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP3);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP3);

        GL.drawElements(GL.TRIANGLES, lefteye_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // MATA KANAN
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP4);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP4);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP4);

        GL.drawElements(GL.TRIANGLES, righteye_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // TELINGA KIRI
        GL.bindBuffer(GL.ARRAY_BUFFER, L_EAR_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, L_EAR_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, L_EAR_FACES);

        GL.drawElements(GL.TRIANGLES, L_ear.faces.length, GL.UNSIGNED_SHORT, 0);

        // TELINGA KANAN
        GL.bindBuffer(GL.ARRAY_BUFFER, R_EAR_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, R_EAR_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, R_EAR_FACES);

        GL.drawElements(GL.TRIANGLES, R_ear.faces.length, GL.UNSIGNED_SHORT, 0);

        // LENGAN KIRI ATAS
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP6);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP6);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP6);

        GL.drawElements(GL.TRIANGLES, upperleft_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // LENGAN KIRI BAWAH
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP7);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP7);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP7);

        GL.drawElements(GL.TRIANGLES, bottomleft_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // LENGAN KANAN ATAS
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP8);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP8);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP8);

        GL.drawElements(GL.TRIANGLES, upperright_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // LENGAN KANAN BAWAH
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP9);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP9);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP9);

        GL.drawElements(GL.TRIANGLES, bottomright_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // KAKI KIRI
        GL.bindBuffer(GL.ARRAY_BUFFER, L_FOOT_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, L_FOOT_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, L_FOOT_FACES);

        GL.drawElements(GL.TRIANGLES, L_foot.faces.length, GL.UNSIGNED_SHORT, 0);

        // KAKI KIRI BAWAH
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP10);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP10);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP10);

        GL.drawElements(GL.TRIANGLES, f_bottom_left_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // KAKI KANAN
        GL.bindBuffer(GL.ARRAY_BUFFER, R_FOOT_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, R_FOOT_COLORS);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, R_FOOT_FACES);

        GL.drawElements(GL.TRIANGLES, R_foot.faces.length, GL.UNSIGNED_SHORT, 0);

        // KAKI KANAN BAWAH
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_bufferP11);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_bufferP11);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_bufferP11);

        GL.drawElements(GL.TRIANGLES, f_bottom_right_indicesP1.length, GL.UNSIGNED_SHORT, 0);

        // ALIS KANAN
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alisKanan_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_alisKanan_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.drawElements(GL.LINE_STRIP, alisKanan.vertices.length, GL.UNSIGNED_SHORT, 0);

        // ALIS KANAN
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alisKiri_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_alisKiri_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.drawElements(GL.LINE_STRIP, alisKiri.vertices.length, GL.UNSIGNED_SHORT, 0);
        
        // HEAD
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_buffer);

        GL.drawElements(GL.TRIANGLES, head_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // JAW
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_jaw_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_jaw_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_jaw_buffer);

        GL.drawElements(GL.TRIANGLES, jaw_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Congor
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_congor_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_congor_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_congor_buffer);

        GL.drawElements(GL.TRIANGLES, congor_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Hidung
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_hidung_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_hidung_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_hidung_buffer);

        GL.drawElements(GL.TRIANGLES, hidung_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Telinga 1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_telinga1_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_telinga1_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_telinga1_buffer);

        GL.drawElements(GL.TRIANGLES, telinga1_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Telinga 1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_telinga2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_telinga2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_telinga2_buffer);

        GL.drawElements(GL.TRIANGLES, telinga2_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mata 1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mata1_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_mata1_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_mata1_buffer);

        GL.drawElements(GL.TRIANGLES, mata1_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mata 2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mata2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_mata2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_mata2_buffer);

        GL.drawElements(GL.TRIANGLES, mata2_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Alis 1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alis1_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_alis1_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_alis1_buffer);

        GL.drawElements(GL.TRIANGLES, alis1_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Alis 2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_alis2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_alis2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_alis2_buffer);

        GL.drawElements(GL.TRIANGLES, alis2_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Badan
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_buffer);

        GL.drawElements(GL.TRIANGLES, badan_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // TAK
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtak_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtak_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtak_buffer);

        GL.drawElements(GL.TRIANGLES, poohtak_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // TAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtakn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtakn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtakn_buffer);

        GL.drawElements(GL.TRIANGLES, poohtakn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // KAK
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkak_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkak_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkak_buffer);

        GL.drawElements(GL.TRIANGLES, poohkak_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // KAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkakn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkakn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkakn_buffer);

        GL.drawElements(GL.TRIANGLES, poohkakn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // tbk
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtbk_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtbk_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtbk_buffer);

        GL.drawElements(GL.TRIANGLES, poohtbk_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // tbkn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohtbkn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohtbkn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohtbkn_buffer);

        GL.drawElements(GL.TRIANGLES, poohtbkn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // kbk
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkbk_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkbk_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkbk_buffer);

        GL.drawElements(GL.TRIANGLES, poohkbk_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // kbkn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_poohkbkn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_poohkbkn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_poohkbkn_buffer);

        GL.drawElements(GL.TRIANGLES, poohkbkn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mulut hehe
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_pooh_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_pooh_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.drawElements(GL.LINE_STRIP, mulut_pooh.vertices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mulut hehe
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_piglet_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_piglet_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.drawElements(GL.LINE_STRIP, mulut_piglet.vertices.length, GL.UNSIGNED_SHORT, 0);

        // RABBIT
        // HEAD
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_head_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_head_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_head_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_head_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // JAW
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_jaw_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_jaw_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_jaw_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_jaw_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Congor
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_congor_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_congor_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_congor_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_congor_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Hidung
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_hidung_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_hidung_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_hidung_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_hidung_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mata luar
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_luar_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_luar_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_luar_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_mata_luar_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // Mata Dalam
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_dalam_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_dalam_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_dalam_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_mata_dalam_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // Mata luar
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_luar_kanan_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_luar_kanan_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_luar_kanan_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_mata_luar_kanan_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // Mata Dalam
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_mata_dalam_kanan_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_mata_dalam_kanan_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_mata_dalam_kanan_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_mata_dalam_kanan_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_telinga
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_telinga_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_telinga_dl
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_dl_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_dl_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_dl_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_telinga_dl_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_telinga_kn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_kn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_kn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_kn_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_telinga_kn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_telinga_kn_dl
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_telinga_kn_dl_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_telinga_kn_dl_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_telinga_kn_dl_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_telinga_kn_dl_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_alis1
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_alis1_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_alis1_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_alis1_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_alis1_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbit_alis2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbit_alis2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbit_alis2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbit_alis2_buffer);

        GL.drawElements(GL.TRIANGLES, rabbit_alis2_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // mulut_rabbit
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_mulut_rabbit_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        
        GL.bindBuffer(GL.ARRAY_BUFFER, color_mulut_rabbit_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        
        GL.drawElements(GL.LINE_STRIP, mulut_rabbit.vertices.length, GL.UNSIGNED_SHORT, 0);
        
        // badan_rabbit
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_buffer);

        GL.drawElements(GL.TRIANGLES, badan_rabbit_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // badan_rabbit2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit2_buffer);

        GL.drawElements(GL.TRIANGLES, badan_rabbit2_indices.length, GL.UNSIGNED_SHORT, 0);

        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, MODEL_MATRIX);
        
        // badan_rabbit_dl
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_dl_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_dl_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_dl_buffer);

        GL.drawElements(GL.TRIANGLES, badan_rabbit_dl_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // badan_rabbit_dl2
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_badan_rabbit_dl2_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_badan_rabbit_dl2_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_badan_rabbit_dl2_buffer);

        GL.drawElements(GL.TRIANGLES, badan_rabbit_dl2_indices.length, GL.UNSIGNED_SHORT, 0);

        // TAK
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittak_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittak_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittak_buffer);

        GL.drawElements(GL.TRIANGLES, rabbittak_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // TAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittakn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittakn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittakn_buffer);

        GL.drawElements(GL.TRIANGLES, rabbittakn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // KAK
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkak_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkak_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkak_buffer);

        GL.drawElements(GL.TRIANGLES, rabbitkak_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // KAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkakn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkakn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkakn_buffer);

        GL.drawElements(GL.TRIANGLES, rabbitkakn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // KAK
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittbk_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittbk_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittbk_buffer);

        GL.drawElements(GL.TRIANGLES, rabbittbk_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // KAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbittbkn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbittbkn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbittbkn_buffer);

        GL.drawElements(GL.TRIANGLES, rabbittbkn_indices.length, GL.UNSIGNED_SHORT, 0);
        
        // rabbitkbkn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkbk_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkbk_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkbk_buffer);

        GL.drawElements(GL.TRIANGLES, rabbitkbk_indices.length, GL.UNSIGNED_SHORT, 0);
       
        // KAKn
        GL.bindBuffer(GL.ARRAY_BUFFER, vertex_rabbitkbkn_buffer);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ARRAY_BUFFER, color_rabbitkbkn_buffer);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);

        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, index_rabbitkbkn_buffer);

        GL.drawElements(GL.TRIANGLES, rabbitkbkn_indices.length, GL.UNSIGNED_SHORT, 0);
        

        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, MODEL_MATRIX);

        // Mendefinisikan titik kontrol untuk kurva Bezier
        
        // Memanggil fungsi untuk setiap kurva

        GL.flush();

        window.requestAnimationFrame(animate);
    };


    animate(0);



}

function mainCanvas(){
    var CANVAS = document.getElementById("mycanvas");
    
    CANVAS.width = window.innerWidth;
    CANVAS.height = window.innerHeight;
    

        // Variabel untuk menyimpan koordinat kotak
        var ctx = CANVAS.getContext("2d");
        var kotak = [];
        var jumlahKlik =  0;
    
        // Event untuk klik mouse
        CANVAS.addEventListener("click", function(event) {
            var x = event.clientX;
            var y = event.clientY;
    
            // Jumlah klik ditambah
            jumlahKlik++;
    
            // Gambar kotak
            ctx.beginPath();
            ctx.rect(x -  5, y -  5,  10,  10);
            ctx.fillStyle = "purple";
            ctx.fill();
            ctx.closePath();
    
            // Simpan koordinat kotak
            kotak.push({ x: x, y: y });
    
            // Gambar garis jika ini klik kedua atau lebih
            if (jumlahKlik >=  2) {
                var kotakSebelumnya = kotak[kotak.length -  2];
                ctx.beginPath();
                ctx.moveTo(kotakSebelumnya.x, kotakSebelumnya.y);
                ctx.lineTo(x, y);
                ctx.strokeStyle = "red";
                ctx.stroke();
                ctx.closePath();
            }
    
            // Gambar kurva jika ini klik ketiga
            if (jumlahKlik >=  3) {
                var kotakPertama = kotak[0];
                var kotakTerakhir = kotak[kotak.length -  1];
    
                ctx.beginPath();
                ctx.moveTo(kotakPertama.x, kotakPertama.y);
                for (var i =  1; i < kotak.length -  1; i++) {
                    // Hitung titik kontrol yang lebih jauh dari garis
                    var titikKontrolX = kotak[i].x + (kotak[i +  1].x - kotak[i].x) /  2;
                    var titikKontrolY = kotak[i].y + (kotak[i +  1].y - kotak[i].y) /  2;
                    ctx.quadraticCurveTo(kotak[i].x, kotak[i].y, titikKontrolX, titikKontrolY);
                }
                ctx.quadraticCurveTo(kotakTerakhir.x, kotakTerakhir.y, kotakTerakhir.x, kotakTerakhir.y);
                ctx.strokeStyle = "grey";
                ctx.stroke();
                ctx.closePath();
            }
        });   
    
}

window.addEventListener('load',main);