uniform vec2 resolution;
uniform float time;
uniform int frame;
uniform float fov;
uniform float raymarchMaximumDistance;
uniform float raymarchPrecision;
uniform vec3 camera;
uniform vec3 target;

uniform samplerCube cubemap;
uniform vec3 anchors[15];

//helpers for view
vec2 squareFrame(vec2 screenSize) {
  vec2 position = 2.0 * (gl_FragCoord.xy / screenSize.xy) - 1.0;
  position.x *= screenSize.x / screenSize.y;
  return position;
}
vec2 squareFrame(vec2 screenSize, vec2 coord) {
  vec2 position = 2.0 * (coord.xy / screenSize.xy) - 1.0;
  position.x *= screenSize.x / screenSize.y;
  return position;
}


mat3 calcLookAtMatrix(vec3 origin, vec3 target, float roll) {
  vec3 rr = vec3(sin(roll), cos(roll), 0.0);
  vec3 ww = normalize(target - origin);
  vec3 uu = normalize(cross(ww, rr));
  vec3 vv = normalize(cross(uu, ww));
  return mat3(uu, vv, ww);
}

//https://github.com/stackgl/glsl-camera-ra
vec3 getRay(mat3 camMat, vec2 screenPos, float lensLength) {
  return normalize(camMat * vec3(screenPos, lensLength));
}

vec3 getRay(vec3 origin, vec3 target, vec2 screenPos, float lensLength) {
  mat3 camMat = calcLookAtMatrix(origin, target, 0.0);
  return getRay(camMat, screenPos, lensLength);
}

mat3 rotationMatrix3(vec3 axis, float angle){
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c          );
}

// http://iquilezles.org/www/articles/smin/smin.htm
float smin( float a, float b, float k ){
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

// http://iquilezles.org/www/articles/smin/smin.htm
vec2 smin( vec2 a, vec2 b, float k ){
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

// http://iquilezles.org/www/articles/smin/smin.htm
float smax( float a, float b, float k ){
    float h = max(k-abs(a-b),0.0);
    return max(a, b) + h*h*0.25/k;
}

vec2 unionAB(vec2 a, vec2 b){return vec2(min(a.x, b.x),1.);}
vec2 intersectionAB(vec2 a, vec2 b){return vec2(max(a.x, b.x),1.);}
vec2 blendAB( vec2 a, vec2 b, float t ){ return vec2(mix(a.x, b.x, t ),1.);}
vec2 subtract(vec2 a, vec2 b){ return vec2(max(-a.x, b.x),1.); }

vec2 sphere( vec3 p, float radius, vec3 pos , vec4 quat){
    mat3 transform = rotationMatrix3( quat.xyz, quat.w );
    float d = length( ( p * transform )-pos ) - radius;
    return vec2(d,0.);
}

vec2 sphere( vec3 p, float radius, vec3 pos ){
    float d = length( p -pos ) - radius;
    return vec2(d,0.);
}

vec2 roundBox(vec3 p, vec3 size, float corner, vec3 pos, vec4 quat ){
    mat3 transform = rotationMatrix3( quat.xyz, quat.w );
    return vec2( length( max( abs( ( p-pos ) * transform )-size, 0.0 ) )-corner,1.);
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdSphere( vec3 p, float s ){
    return length(p)-s;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdEllipsoid( vec3 p, vec3 r ) {// approximated
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

vec2 sdStick(vec3 p, vec3 a, vec3 b, float r1, float r2){ // approximated
    vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return vec2( length( pa - ba*h ) - mix(r1,r2,h*h*(3.0-2.0*h)), h );
}

vec2 line( vec3 p, vec3 a, vec3 b, float r ){
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return vec2( length( pa - ba*h ) - r, 1. );
}

vec2 plane( vec3 p , vec3 n) { return vec2( dot(p, n), 1. ); }
//with height
vec2 plane( vec3 p , vec4 n) { return vec2( dot(p, n.xyz) + n.w, 1. ); }

// http://iquilezles.org/www/articles/smin/smin.htm
vec4 opU( vec4 d1, vec4 d2 ){
	return (d1.x<d2.x) ? d1 : d2;
}


//http://www.pouet.net/topic.php?post=367360
const vec3 pa = vec3(1., 57., 21.);
const vec4 pb = vec4(0., 57., 21., 78.);
float perlin(vec3 p) {
	vec3 i = floor(p);
	vec4 a = dot( i, pa ) + pb;
	vec3 f = cos((p-i)*acos(-1.))*(-.5)+.5;
	a = mix(sin(cos(a)*a),sin(cos(1.+a)*(1.+a)), f.x);
	a.xy = mix(a.xz, a.yw, f.y);
	return mix(a.x, a.y, f.z);
}

float zigzag( float x, float m ){
    return abs( mod( x, (2.*m) ) -m);
}

float href;
float hsha;

vec4 map( vec3 position){
vec4 quat = vec4( 1., 0., 0., .5 );
    float rad = 500.;
    vec3 dir = vec3(.0,.0, time * 4.);
    vec2 ground = sphere( position + perlin( ( position + dir ) * .1 ), rad, vec3( 0.,-rad + 2.,0. ) );
    //ground = unionAB( ground, plane( position - vec3( 0.,100.,0. ), vec3( 0.,-1.,0. ) ) );

    float o = zigzag( position.x, .25 ) + zigzag( position.x, .21 );

    float radius = .35;
    float blendFactor = 1.;
    dir = vec3( 0., -time * 3., 0. );

    float s = fract( sin( sin( floor( position.x / 0.01 ) * 2. ) / 0.01 ) * 10. ) * 0.;

    vec2 skeleton = line( position, anchors[0] + vec3(0.,1.,0.), anchors[1], .5 ) + perlin(position+dir )*1.5;

    skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*1.2 , 1. ) );

    //blend distance (color blend)
    float dis0 = skeleton.x;

    //left arm
    skeleton = smin( skeleton, line( position, anchors[1], anchors[2], radius ), blendFactor );//shoulder L
    skeleton = smin( skeleton, line( position, anchors[2], anchors[3], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[3], anchors[4], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[4], quat ) + o, blendFactor );

    //right arm
    skeleton = smin( skeleton, line( position, anchors[1], anchors[5], radius ), blendFactor );//shoulder R
    skeleton = smin( skeleton, line( position, anchors[5], anchors[6], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[6], anchors[7], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[7], quat ) + o, blendFactor );

    //spine
    skeleton = smin( skeleton, line( position, anchors[1], anchors[8], radius * 1.5 ), blendFactor );

    //belly
    skeleton = smin( skeleton, sphere( position, radius * 2., anchors[8] ), blendFactor );

    //left leg
    skeleton = smin( skeleton, line( position, anchors[8], anchors[9], radius * .5 ), blendFactor * .75 );

    skeleton = smin( skeleton, line( position, anchors[9], anchors[10], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[10], anchors[11], radius ), blendFactor );

    //right leg
    skeleton = smin( skeleton, line( position, anchors[8], anchors[12], radius * .5 ), blendFactor * .75 );

    skeleton = smin( skeleton, line( position, anchors[12], anchors[13], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[13], anchors[14], radius ), blendFactor * 1.5 );

    vec2 _out = smin( ground, skeleton, 1. );
    _out.y = smoothstep( 0., dis0, _out.x );
    
    return vec4(_out, 1.0, 1.0);
}

vec4 castRay( vec3 ro, vec3 rd){
    vec4 res = vec4(-1.0,-1.0,0.0,1.0);

    float tmin = 0.0;
    float tmax = raymarchMaximumDistance;
    
    float t = tmin;
    for( int i=0; i<256; i++ ){
    	if(t<tmax){
    		vec4 h = map( ro+rd*t );
	        if( abs(h.x)<(raymarchPrecision) ){ 
	            res = h;//vec4(t,h.yzw); 
	            break;
	        }
	        t += h.x;
    	}
    }
    
    return res;
}


float calcSoftshadow(vec3 ro, vec3 rd, float time ){
    float res = 1.0;
    float tmax = 12.0;
    
    float t = 0.02;
    for( int i=0; i<50; i++ ){
		float h = map( ro + rd*t).x;
        res = min( res, mix(1.0,16.0*h/t, hsha) );
        t += clamp( h, 0.05, 0.40 );
        if( res<0.005 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

vec3 calcNormal( vec3 pos  ){

    const vec3 v1 = vec3( 1.0,-1.0,-1.0);
  const vec3 v2 = vec3(-1.0,-1.0, 1.0);
  const vec3 v3 = vec3(-1.0, 1.0,-1.0);
  const vec3 v4 = vec3( 1.0, 1.0, 1.0);

  return normalize( v1 * map( pos + v1*0.002 ).x +
                    v2 * map( pos + v2*0.002 ).x +
                    v3 * map( pos + v3*0.002 ).x +
                    v4 * map( pos + v4*0.002 ).x );
}

float calcOcclusion( vec3 pos, vec3 nor ){
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ ){
        float h = 0.01 + 0.11*float(i)/4.0;
        vec3 opos = pos + h*nor;
        float d = map( opos ).x;
        occ += (h-d)*sca;
        sca *= 0.95;
    }

    return clamp( 1.0 - 2.0*occ, 0.0, 1.0 );
}

vec3 render( vec3 ro,vec3 rd ){
	vec3 col = vec3(0.5, 0.2, 0.9) - max(rd.y,0.0)*0.5;
    // sky clouds
    vec2 uv = 1.5*rd.xz/rd.y;
    float cl  = 1.0*(sin(uv.x)+sin(uv.y)); uv *= mat2(0.8,0.6,-0.6,0.8)*2.1;
          cl += 0.5*(sin(uv.x)+sin(uv.y));
    col += 0.1*(-1.0+2.0*smoothstep(-0.1,0.1,cl-0.4));
    // sky horizon
	col = mix( col, vec3(0.5, 0.7, .9), exp(-10.0*max(rd.y,0.0)) );    

	vec4 res = castRay(ro,rd);

	return res.xyz;
} 

void main(){
	vec2  screenPos    = squareFrame( resolution );
    float alpha = smoothstep( .2, 1., 1.- abs( screenPos.x ) );
    //uncomment below for fullscreen
    alpha = 1.;
    if( alpha <= 0. ) discard;

    vec3 rd = getRay( camera, target, screenPos, fov );
    vec3 col = render( camera, rd);
    /*vec3 co = vec3(0.8, 0.3, 0.3) - max(rd.y,0.0)*0.5;

    gl_FragColor = vec4(mix( col, vec3(1.), screenPos.y), alpha );
*/
    // color grading
    /*col = col*vec3(1.11,0.89,0.79);
    // compress        
    col = 1.35*col/(1.0+col);
    // gamma
    col = pow( col, vec3(0.4545) );*/

    gl_FragColor = vec4(col, 1.0);
}