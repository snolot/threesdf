
uniform vec2 resolution;
uniform float time;
uniform float fov;
uniform float raymarchMaximumDistance;
uniform float raymarchPrecision;
uniform vec3 camera;
uniform vec3 target;

uniform samplerCube cubemap;
uniform vec3 anchors[15];


#define MOD3 vec3(.16532,.17369,.15787)
#define MOD2 vec2(.16632,.17369)

float hash12(vec2 p)
{
    p  = fract(p * MOD2);
    p += dot(p.xy, p.yx+19.19);
    return fract(p.x * p.y);
}

float Hash(vec3 p){
    p  = fract(p * MOD3);
    p += dot(p.xyz, p.yzx + 19.19);
    return fract(p.x * p.y * p.z);
}

float Noise3d(vec3 p)
{
    vec2 add = vec2(1.0, 0.0);
    p *= 10.0;
    float h = 0.0;
    float a = .3;
    for (int n = 0; n < 4; n++)
    {
        vec3 i = floor(p);
        vec3 f = fract(p); 
        f *= f * (3.0-2.0*f);

        h += mix(
            mix(mix(Hash(i), Hash(i + add.xyy),f.x),
                mix(Hash(i + add.yxy), Hash(i + add.xxy),f.x),
                f.y),
            mix(mix(Hash(i + add.yyx), Hash(i + add.xyx),f.x),
                mix(Hash(i + add.yxx), Hash(i + add.xxx),f.x),
                f.y),
            f.z)*a;
         a*=.5;
        p += p;
    }
    return h;
}
//uses most of the StackGL methods
//https://github.com/stackgl

//https://github.com/hughsk/glsl-square-frame

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

//https://github.com/stackgl/glsl-look-at/blob/gh-pages/index.glsl

mat3 calcLookAtMatrix(vec3 origin, vec3 target, float roll) {
  vec3 rr = vec3(sin(roll), cos(roll), 0.0);
  vec3 ww = normalize(target - origin);
  vec3 uu = normalize(cross(ww, rr));
  vec3 vv = normalize(cross(uu, ww));
  return mat3(uu, vv, ww);
}

//https://github.com/stackgl/glsl-camera-ray

vec3 getRay(mat3 camMat, vec2 screenPos, float lensLength) {
  return normalize(camMat * vec3(screenPos, lensLength));
}
vec3 getRay(vec3 origin, vec3 target, vec2 screenPos, float lensLength) {
  mat3 camMat = calcLookAtMatrix(origin, target, 0.0);
  return getRay(camMat, screenPos, lensLength);
}

/////////////////////////////////////////////////////////////////////////

mat3 rotationMatrix3(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c          );
}

/////////////////////////////////////////////////////////////////////////

//primitives

vec2 sphere( vec3 p, float radius, vec3 pos , vec4 quat)
{
    mat3 transform = rotationMatrix3( quat.xyz, quat.w );
    float d = length( ( p * transform )-pos ) - radius;
    return vec2(d,0.);
}

vec2 sphere( vec3 p, float radius, vec3 pos )
{
    float d = length( p -pos ) - radius;
    return vec2(d,0.);
}

vec2 roundBox(vec3 p, vec3 size, float corner, vec3 pos, vec4 quat )
{
    mat3 transform = rotationMatrix3( quat.xyz, quat.w );
    return vec2( length( max( abs( ( p-pos ) * transform )-size, 0.0 ) )-corner,1.);
}

vec2 line( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return vec2( length( pa - ba*h ) - r, 1. );
}
//operations

vec2 unionAB(vec2 a, vec2 b){return vec2(min(a.x, b.x),1.);}
vec2 intersectionAB(vec2 a, vec2 b){return vec2(max(a.x, b.x),1.);}
vec2 blendAB( vec2 a, vec2 b, float t ){ return vec2(mix(a.x, b.x, t ),1.);}
vec2 subtract(vec2 a, vec2 b){ return vec2(max(-a.x, b.x),1.); }
//http://iquilezles.org/www/articles/smin/smin.htm
vec2 smin( vec2 a, vec2 b, float k ) { float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 ); return vec2( mix( b.x, a.x, h ) - k*h*(1.0-h), 1. ); }
float smin( float a, float b, float k ) { float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 ); return mix( b, a, h ) - k*h*(1.0-h); }

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
float zigzag( float x, float m )
{
    return abs( mod( x, (2.*m) ) -m);
}
/////////////////////////////////////////////////////////////////////////

// STOP ! ! !

// HAMMER TIME !

/////////////////////////////////////////////////////////////////////////
//AshimaOptim https://www.shadertoy.com/view/Xd3GRf
vec4 permute(vec4 x){return mod(x*x*34.0+x,289.);}
float snoise(vec3 v){
  const vec2  C = vec2(0.166666667, 0.33333333333) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 = v - i + dot(i, C.xxx) ;
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy;
  vec3 x3 = x0 - D.yyy;
  i = mod(i,289.);
  vec4 p = permute( permute( permute(
      i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
    + i.y + vec4(0.0, i1.y, i2.y, 1.0 ))
    + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));
  vec3 ns = 0.142857142857 * D.wyz - D.xzx;
  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);
  vec4 x_ = floor(j * ns.z);
  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = floor(j - 7.0 * x_ ) *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);
  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));
  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;
  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);
  vec4 norm = inversesqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m * m;
  return .5 + 12.0 * dot( m, vec4( dot(p0,x0), dot(p1,x1),dot(p2,x2), dot(p3,x3) ) );
}

float fogmap(in vec3 p, in float d)
{
    p.xz -= time*2.+sin(p.z*.3)*3.;
    p.y -= time*.3;

    return (max(Noise3d(p*.1+ 1.)-.1,0.0)*Noise3d(p*.02))*.6;
}

const int raymarchSteps = 50;
const float PI = 3.14159;

//no height
vec2 plane( vec3 p , vec3 n) { return vec2( dot(p, n), 1. ); }
//with height
vec2 plane( vec3 p , vec4 n) { return vec2( dot(p, n.xyz) + n.w, 1. ); }
vec2 field( vec3 position ){
    //position
    //vec3 zero = vec3(0.);
    if(position.y < 2.5)
    //position += vec3(fogmap(position, position.x) * 3.);
    position -= snoise(position*0.85) ;
    //rotation
    vec4 quat = vec4( 1., 0., 0., .5 );
    float rad = 100.;
    vec3 dir = vec3(.0,.0, time * 2.);
    vec2 ground = sphere( position + perlin( ( position + dir ) * .1 ), rad, vec3( 0.,-rad + 2.,0. ) );
    ground = unionAB( ground, plane( position - vec3( 0.,100.,0. ), vec3( 0.,-1.,0. ) ) );

    float o = zigzag( position.x, .05 );// + zigzag( position.x, .01 );

    float radius = .75;
    float blendFactor = .8;
    dir = vec3( 0., -time * 1., 0. );

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
    skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );

    //right arm
    skeleton = smin( skeleton, line( position, anchors[1], anchors[5], radius ), blendFactor );//shoulder R
    skeleton = smin( skeleton, line( position, anchors[5], anchors[6], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[6], anchors[7], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[7], quat ) + o, blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*1.2 , 1. ) );
    //spine
    skeleton = smin( skeleton, line( position, anchors[1], anchors[8], radius * 1.5 ), blendFactor );

    //belly
    skeleton = smin( skeleton, sphere( position, radius * 2., anchors[8] ), blendFactor );
    skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );
    //left leg
    skeleton = smin( skeleton, line( position, anchors[8], anchors[9], radius * .5 ), blendFactor * .75 );

    skeleton = smin( skeleton, line( position, anchors[9], anchors[10], radius ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.8 , 1. ) );
    skeleton = smin( skeleton, line( position, anchors[10], anchors[11], radius ), blendFactor );

    //right leg
    skeleton = smin( skeleton, line( position, anchors[8], anchors[12], radius * .5 ), blendFactor * .75 );

    skeleton = smin( skeleton, line( position, anchors[12], anchors[13], radius ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );
    skeleton = smin( skeleton, line( position, anchors[13], anchors[14], radius ), blendFactor * 1.5 );

    vec2 _out = smin( ground, skeleton, 1. );
    _out.y = smoothstep( 0., dis0, _out.x );
    
    return _out;
}

/////////////////////////////////////////////////////////////////////////

// the methods below this need the field function

/////////////////////////////////////////////////////////////////////////


//the actual raymarching from:
//https://github.com/stackgl/glsl-raytrace/blob/master/index.glsl
vec2 raymarching( vec3 rayOrigin, vec3 rayDir, float maxd, float precis ) {

    float latest = precis;
    float dist   = 0.0;
    float type   = -1.0;
    for (int i = 0; i < raymarchSteps; i++) {

        if (latest < precis || dist > maxd) break;

        vec2 result = field( rayOrigin + rayDir * dist );

        latest = result.x;
        dist  += latest;
        type = result.y;
    }

    vec2 res    = vec2(-1.0, -1.0 );
    if (dist < maxd) { res = vec2( dist, type ); }
    return res;

}

//https://github.com/stackgl/glsl-sdf-normal
vec3 calcNormal(vec3 pos, float eps) {
  const vec3 v1 = vec3( 1.0,-1.0,-1.0);
  const vec3 v2 = vec3(-1.0,-1.0, 1.0);
  const vec3 v3 = vec3(-1.0, 1.0,-1.0);
  const vec3 v4 = vec3( 1.0, 1.0, 1.0);

  return normalize( v1 * field( pos + v1*eps ).x +
                    v2 * field( pos + v2*eps ).x +
                    v3 * field( pos + v3*eps ).x +
                    v4 * field( pos + v4*eps ).x );
}

vec3 calcNormal(vec3 pos) {
  return calcNormal(pos, 0.002);
}

vec3 rimlight( vec3 pos, vec3 nor ){
    vec3 v = normalize(-pos);
    float vdn = 1.0 - max(dot(v, nor), 0.0);
    return vec3(smoothstep(0., 1.0, vdn));
}

vec3 phongColor(in vec3 p, in vec3 rd){
    
    vec3 Ks = vec3(1.0); // specular reflected intensity
    float shininess = 80.0;
    
    vec3 n = calcNormal(p, 0.01);
    vec3 ref = reflect( rd, n );
    vec3 final = vec3(0);
    

    // light 0
    {
        vec3 light_pos   = p+vec3( -1.0, .3, -1.0 );
        vec3 light_color = vec3( 1.0, 1.0, 1.0 );
    
        vec3 vl = normalize( light_pos - p );
        vec3 specular = vec3( max( 0.0, dot( vl, ref ) ) );
        specular = pow( specular, vec3( shininess ) );
        final += light_color * specular; 
    }
    return final;
}


/////////////////////////////////////////////////////////////////////////
// SHADING
/////////////////////////////////////////////////////////////////////////

float checkeredPattern( vec3 p ) {
    float u = 1.0 - floor( mod( p.x, 2.0 ) );
    float v = 1.0 - floor( mod( p.z, 2.0 ) );

    if ( ( u == 1.0 && v < 1.0 ) || ( u < 1.0 && v == 1.0 ) ) {
        return 0.2;
    } else {
        return 1.0;
    }
}

vec3 hsv2rgb( vec3 c ) {
    vec4 K = vec4( 1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0 );
    vec3 p = abs( fract( c.xxx + K.xyz ) * 6.0 - K.www );
    return c.z * mix( K.xxx, clamp( p - K.xxx, 0.0, 1.0 ), c.y );
}


void main() {
    vec2  screenPos    = squareFrame( resolution );
    float alpha = smoothstep( .2, 1., 1.- abs( screenPos.x ) );
    //uncomment below for fullscreen
    alpha = 1.;
    if( alpha <= 0. ) discard;

    vec3  rayDirection = getRay( camera, target, screenPos, fov );
    vec2 collision = raymarching( camera, rayDirection, raymarchMaximumDistance, raymarchPrecision );
    
    //vec3 col = vec3( 0.85 );
    vec3 col = vec3(0.5, 0.6, 0.8) - max(rayDirection.y-.3,0.0)*0.45;
    // sky clouds
    /*vec2 uv = 1.5*rayDirection.xz/rayDirection.y;
    float cl  = 1.0*(sin(uv.x)+sin(uv.y));
    uv *= mat2(0.8,0.6,-0.6,0.8)*2.1;
    cl += 0.5*(sin(uv.x)+sin(uv.y));
    col += 0.1*(-1.0+2.0*smoothstep(-0.1,0.1,cl-0.4));
    // sky horizon*/
    //col = mix( col, vec3(0.5, 0.7, .9), exp(-10.0*max(rayDirection.y,0.0)) );    
    

    //gl_FragColor = vec4(mix( col, vec3(1.), screenPos.y), alpha );

    
    if ( collision.x > -0.5){

        vec3 pos = camera + rayDirection * collision.x;
        vec3 nor = calcNormal( pos,.1 );
        vec3 tex = textureCube( cubemap, nor ).rgb * vec3(1.8, 1.6, 0.7);
        /*col = mix( col, tex, collision.y );
        col = col * rimlight( pos, nor )  + tex * .5;*/
        //col = col + tex * .2;
        col = mix(phongColor(pos, nor),tex, .9) ;
        gl_FragColor = vec4( col * rimlight( pos, nor ), alpha );
    }
}