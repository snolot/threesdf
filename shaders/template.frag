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

float hash12(vec2 p){
    p  = fract(p * MOD2);
    p += dot(p.xy, p.yx+19.19);
    return fract(p.x * p.y);
}

float Hash(vec3 p){
    p  = fract(p * MOD3);
    p += dot(p.xyz, p.yzx + 19.19);
    return fract(p.x * p.y * p.z);
}

float Noise3d(in vec3 p){
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

mat3 rotationMatrix(vec3 axis, float angle){
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

mat3 calcLookAtMatrix(vec3 origin, vec3 target, float roll) {
  vec3 rr = vec3(sin(roll), cos(roll), 0.0);
  vec3 ww = normalize(target - origin);
  vec3 uu = normalize(cross(ww, rr));
  vec3 vv = normalize(cross(uu, ww));
  return mat3(uu, vv, ww);
}

vec3 getRay(mat3 camMat, vec2 screenPos, float lensLength) {
  return normalize(camMat * vec3(screenPos, lensLength));
}

vec3 getRay(vec3 origin, vec3 target, vec2 screenPos, float lensLength) {
  mat3 camMat = calcLookAtMatrix(origin, target, 0.0);
  return getRay(camMat, screenPos, lensLength);
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

vec2 sphere( vec3 p, float radius, vec3 pos , vec4 quat){
    mat3 transform = rotationMatrix( quat.xyz, quat.w );
    float d = length( ( p * transform )-pos ) - radius;
    return vec2(d,0.);
}

vec2 sphere( vec3 p, float radius, vec3 pos ){
    float d = length( p -pos ) - radius;
    return vec2(d,0.);
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
vec2 sdSphere( vec3 p, float s, vec3 pos){
    return vec2(length(p-pos)-s, 0.);
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdEllipsoid( in vec3 p, in vec3 r ){ // approximated
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdTorus( vec3 p, vec2 t ){
  return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

vec2 sdCapsule( vec3 p, vec3 a, vec3 b, float r ){
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return vec2(length( pa - ba*h ) - r, 1.);
	
}

vec2 sdStick(vec3 p, vec3 a, vec3 b, float r1, float r2){ // approximated
    vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return vec2( length( pa - ba*h ) - mix(r1,r2,h*h*(3.0-2.0*h)), h );
}

vec2 roundBox(vec3 p, vec3 size, float corner, vec3 pos, vec4 quat ){
    mat3 transform = rotationMatrix( quat.xyz, quat.w );
    return vec2( length( max( abs( ( p-pos ) * transform )-size, 0.0 ) )-corner,1.);
}

float lengthN(in vec2 p, in float n)
{
	p = pow(abs(p), vec2(n));
	return pow(p.x+p.y, 1.0/n);
}

vec4 opU( vec4 d1, vec4 d2 ){
	return (d1.x<d2.x) ? d1 : d2;
}

vec2 blendAB( vec2 a, vec2 b, float t ){ 
	return vec2(mix(a.x, b.x, t ),1.);
}

vec2 intersectionAB(vec2 a, vec2 b){
	return vec2(max(a.x, b.x),1.);
}

vec2 unionAB(vec2 a, vec2 b){
	return vec2(min(a.x, b.x),1.);
}


vec2 field( vec3 position ){
    
    vec4 quat = vec4( 1., 0., 0., .5 );
    float rad = 100.;
    vec3 dir = vec3(.0,time * .1, .0);
    
    float radius = .55;
    float blendFactor = 0.0;
    
    //vec2 skeleton = sdCapsule( position, anchors[0] + vec3(0.,1.,0.), anchors[1], .8 ) ;
    vec2 skeleton = sdSphere( position, 1.2 , anchors[0] + vec3(0.,1.,0.)) ;
    //skeleton = blendAB( skeleton, vec2( perlin(position+dir)*.2 , .1 ) , .65 );

    skeleton.y = 2.;

    float dis0 = skeleton.x;

    //left arm
    //skeleton = smin( skeleton, sdCapsule( position, anchors[1], anchors[2], radius ), blendFactor );//shoulder L
    skeleton = smin( skeleton, sdCapsule( position, anchors[2], anchors[3], radius ), blendFactor );
    skeleton = smin( skeleton, sdCapsule( position, anchors[3], anchors[4], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[4], quat ), blendFactor );
    //skeleton = blendAB( skeleton, vec2( perlin(position+dir)*.2 , .1 ), .55 );

    //right arm
    skeleton = smin( skeleton, sdCapsule( position, anchors[2], anchors[5], radius ), blendFactor );//shoulder R
    //skeleton = smin( skeleton, sdCapsule( position, anchors[1], anchors[5], radius ), blendFactor );//shoulder R
    skeleton = smin( skeleton, sdCapsule( position, anchors[5], anchors[6], radius ), blendFactor );
    skeleton = smin( skeleton, sdCapsule( position, anchors[6], anchors[7], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[7], quat ), blendFactor );
    //skeleton = blendAB( skeleton, vec2( perlin(position+dir)*.2 , .1 ) , .55 );
    
    //spine
    //skeleton = smin( skeleton, sdCapsule( position, anchors[8], (anchors[9] + ((anchors[12] - anchors[9]) * vec3(.5))), radius * 1. ), blendFactor );
    skeleton = smin( skeleton, sdCapsule( position, (anchors[2] + ((anchors[5] - anchors[2]) * vec3(.5))), anchors[8], radius * 1. ), blendFactor );

    //belly
    //skeleton = smin( skeleton, sdCapsule( position, anchors[8], anchors[12] , radius * 1. ), blendFactor );
    //skeleton = smin( skeleton, sphere( position, radius * 2., anchors[8] ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );
    //left leg
    skeleton = smin( skeleton, sdCapsule( position, anchors[8], (anchors[9] + ((anchors[12] - anchors[9]) * vec3(.5))), radius * 1. ), blendFactor );

    skeleton = smin( skeleton, sdCapsule( position, anchors[9], anchors[10], radius ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.8 , 1. ) );
    skeleton = smin( skeleton, sdCapsule( position, anchors[10], anchors[11], radius ), blendFactor );

    //right leg
    skeleton = smin( skeleton, sdCapsule( position, anchors[9], anchors[12], radius * 1. ), blendFactor);

    skeleton = smin( skeleton, sdCapsule( position, anchors[12], anchors[13], radius ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );
    skeleton = smin( skeleton, sdCapsule( position, anchors[13], anchors[14], radius ), blendFactor * 1. );
    //skeleton = blendAB( skeleton, vec2( perlin(position+dir)*.05 , .01 ) , .65 );
    //vec2 _out = smin( ground, skeleton, 1. );
    //skeleton.y = smoothstep( 0., dis0, skeleton.x );
    

    return skeleton;
}

const int raymarchSteps = 50;

vec2 raymarching( vec3 rayOrigin, vec3 rayDir, float maxd, float precis ) {

    float latest = precis * 2.0;
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

float hash( float n ){
	return fract(sin(n)*3538.5453);
}

float calcSoftshadow( in vec3 ro, in vec3 rd ) {
  float res = 1.0;
  float t = 0.0005;                 // selfintersection avoidance distance
  float h = 1.0;
  for( int i=0; i<40; i++ ) {
    h = field(ro + rd*t).x;
    res = min( res, 64.0*h/t );   // 64 is the hardness of the shadows
    t += clamp( h, 0.02, 2.0 );   // limit the max and min stepping distances
  }
  return clamp(res,0.0,1.0);
}

float calcAO( in vec3 p, in vec3 n, float maxDist, float falloff ){
	float ao = 0.0;
	const int nbIte = 6;
	for( int i=0; i<nbIte; i++ )
	{
		float l = hash(float(i))*maxDist;
		vec3 rd = n*l;
		ao += (l - field( p + rd ).x) / pow(1.+l, falloff);
	}
	return clamp( 1.-ao/float(nbIte), 0., 1.);
}

float thickness( in vec3 p, in vec3 n, float maxDist, float falloff )
{
	float ao = 0.0;
	const int nbIte = 6;
	for( int i=0; i<nbIte; i++ )
	{
		float l = hash(float(i))*maxDist;
		vec3 rd = -n*l;
		ao += (l + field( p + rd ).x) / pow(1.+l, falloff);
	}
	return clamp( 1.-ao/float(nbIte), 0., 1.);
}

vec3 postEffects( in vec3 col, in vec2 uv, in float time )
{
	// gamma correction
	//col = pow( clamp(col,0.0,1.0), vec3(0.45) );
	// vigneting
	//col *= 0.7+.8*pow( 16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y), 1.95 );
	col*=1.-pow(length(uv*uv*uv*uv)*1.1,6.);
	return col;
}

void main() {
    vec2  screenPos    = squareFrame( resolution );
    float alpha = smoothstep( .2, 1., 1.- abs( screenPos.x ) );
    
    alpha = 1.;
    if( alpha <= 0.)
    	discard;

    vec3  rayDirection = getRay( camera, target, screenPos, fov );
    vec2 collision = raymarching( camera, rayDirection, raymarchMaximumDistance, raymarchPrecision );
    
    
    vec3 col = vec3(0.7, 0.6, 0.8);// - max(rayDirection.y-.3,0.0)*0.45;
   	
   	vec3 pos = camera + rayDirection * collision.x;
    vec3 nor = calcNormal( pos,.1 );
    vec3 tex = textureCube( cubemap, nor ).rgb;// * vec3(2.0, 1.4, 0.);
   		
    if ( collision.x > -0.5){

   		float ao = calcAO(pos,nor,30.,.9);
   		float thi = thickness(pos, nor, 1., 1.5);
   		
   		vec3 lpos1 = anchors[0] + vec3(0.,1.,0.);//vec3(0.0,35.+sin(time * 2.)*20.,0.0);
		vec3 ldir1 = normalize(lpos1-pos);
		float latt1 = pow( length(lpos1-pos)*.1, 1.5 );
	    float trans1 =  pow( clamp( dot(-rayDirection, -ldir1+nor), 0., 1.), 1.) + 1.;
		vec3 diff1 = vec3(.0,.5,1.) * (max(dot(nor,ldir1),0.) ) / latt1;
		col *=  diff1;
		col += vec3(.3,.2,.05) * (trans1/latt1)*thi;
		col *= ao;
		col = max(vec3(.05),col);

		#define POW2(a) (a*a)

		#define CAPSULECOLOR vec3(1,0,.3)//(vec3(1,0,.3)*(-cos(time*.2)*0.5+0.5))
		col += ( 1.0/POW2(distance(vec3(0.), pos)) ) * 25.0 * CAPSULECOLOR;
    	//col += ( 1.0/POW2(distance(vec3(0.), pos)) ) * 10.0 * CAPSULECOLOR;

    	#define RIMCOLOR vec3(0,0.1,0.3) * max(0.0, sin(atan(float(pos.z), float(pos.x))*5.0+time*5.0)) * step(10.0, length(pos)) * (1.0-smoothstep(1., 0., abs(float(pos.y))))
    	col += clamp(( 1.0/abs(sdTorus(vec3(pos - vec3(0,0,0)), vec2(1.,20.0)) )), 0., 1.0) * 5.0 * RIMCOLOR;
        
    	//col = postEffects( col, screenPos, time );
        gl_FragColor = vec4( col * rimlight( pos, nor ) * 1.0, 1. );
    }else{
    	#define RIMCOLORBG vec3(0,0.1,0.3) * max(0.0, sin(atan(float(pos.z), float(pos.x))*5.0+time*5.0)) * step(1.0, length(pos)) * (1.0-smoothstep(1., 0., abs(float(pos.y))))
    	col += clamp(( 1.0/abs(sdTorus(vec3(pos - vec3(0,0,0)), vec2(10.,24.0)) )), 0., 1.0) * 10.0 * RIMCOLORBG;
        gl_FragColor = vec4( col /** rimlight( pos, nor ) * tex*/, 1. );
    	//discard;
    	/*col = mix( col, vec3(0.0, 0.0, 0.0), exp(-10.0*max(rayDirection.y,0.0)) ); 
    	vec2 uv = 1.5*rayDirection.xz/rayDirection.y;
	    float cl  = 1.0*(sin(uv.x)+sin(uv.y));
	    uv *= mat2(0.8,0.6,-0.6,0.8)*2.1;
	    cl += 0.5*(sin(uv.x)+sin(uv.y));
	    col += 0.1*(-1.0+2.0*smoothstep(-0.1,0.1,cl-0.4));
	    gl_FragColor = vec4( col, 1.); */
    }

}
