
uniform vec2 resolution;
uniform float time;
uniform float fov;
uniform float raymarchMaximumDistance;
uniform float raymarchPrecision;
uniform vec3 camera;
uniform vec3 target;

uniform samplerCube cubemap;
uniform vec3 anchors[15];
// Dither the entire screen for a fun effect
//#define DITHERING
// Whether you want 
//#define CAMERAROTATING
mat3 rotationMatrix(vec3 axis, float angle){
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
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

const int MAX_RAY_STEPS = 90;

float sdSphere(vec3 p, float d) { return length(p) - d; } 

float sdBox( vec3 p, vec3 b ){
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) +length(max(d,0.0));
}

float sdTorus( vec3 p, vec2 t ){
  return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r ){
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
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

vec2 roundBox(vec3 p, vec3 size, float corner, vec3 pos, vec4 quat ){
    mat3 transform = rotationMatrix( quat.xyz, quat.w );
    return vec2( length( max( abs( ( p-pos ) * transform )-size, 0.0 ) )-corner,1.);
}

vec2 line( vec3 p, vec3 a, vec3 b, float r ){
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return vec2( length( pa - ba*h ) - r, 1. );
}

vec2 smin( vec2 a, vec2 b, float k ) { float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 ); return vec2( mix( b.x, a.x, h ) - k*h*(1.0-h), 1. ); }
float smin( float a, float b, float k ) { float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 ); return mix( b, a, h ) - k*h*(1.0-h); }

const int raymarchSteps = 50;
const float PI = 3.14159;

//no height
vec2 plane( vec3 p , vec3 n) { return vec2( dot(p, n), 1. ); }
//with height
vec2 plane( vec3 p , vec4 n) { return vec2( dot(p, n.xyz) + n.w, 1. ); }
vec2 field( vec3 position ){
    //position
    //vec3 zero = vec3(0.);
    //if(position.y < 2.5)
    //position += vec3(fogmap(position, position.x) * 3.);
    //position -= snoise(position*0.85) ;
    //rotation
    vec4 quat = vec4( 1., 0., 0., .5 );
    /*
    float rad = 100.;
    vec3 dir = vec3(.0,.0, time * 2.);
    vec2 ground = sphere( position + perlin( ( position + dir ) * .1 ), rad, vec3( 0.,-rad + 2.,0. ) );
    ground = unionAB( ground, plane( position - vec3( 0.,100.,0. ), vec3( 0.,-1.,0. ) ) );*/

    //float o = zigzag( position.x, .05 );// + zigzag( position.x, .01 );

    float radius = .7;
    float blendFactor = .8;
    //dir = vec3( 0., -time * 1., 0. );

    float s = fract( sin( sin( floor( position.x / 0.01 ) * 2. ) / 0.01 ) * 10. ) * 0.;

   
    vec2 skeleton = line( position, anchors[0] + vec3(0.,1.,0.), anchors[1], .5 );// + perlin(position+dir )*1.5;

    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*1.2 , 1. ) );

    //blend distance (color blend)
    float dis0 = skeleton.x;

    //left arm
    skeleton = smin( skeleton, line( position, anchors[1], anchors[2], radius ), blendFactor );//shoulder L
    skeleton = smin( skeleton, line( position, anchors[2], anchors[3], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[3], anchors[4], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[4], quat ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );

    //right arm
    skeleton = smin( skeleton, line( position, anchors[1], anchors[5], radius ), blendFactor );//shoulder R
    skeleton = smin( skeleton, line( position, anchors[5], anchors[6], radius ), blendFactor );
    skeleton = smin( skeleton, line( position, anchors[6], anchors[7], radius ), blendFactor );

    //hand
    skeleton = smin( skeleton, roundBox( position, vec3( .1,.5,.1 ), .5, anchors[7], quat ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*1.2 , 1. ) );
    //spine
    skeleton = smin( skeleton, line( position, anchors[1], anchors[8], radius * 1.5 ), blendFactor );

    //belly
    skeleton = smin( skeleton, sphere( position, radius * 2., anchors[8] ), blendFactor );
    //skeleton = intersectionAB( skeleton, vec2( perlin(position+dir)*.2 , 1. ) );
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

    //vec2 _out = smin( ground, skeleton, 1. );
    //skeleton.y = smoothstep( 0., dis0, skeleton.x );
    
    return skeleton;
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


// this function returns true if there's a bool in the provided grid position.
bool getVoxel(ivec3 c, mat3 rotMat1, mat3 rotMat2) {
	vec3 p = vec3(c);// + vec3(0.5);
    // Generate 2 rotation matrices for the 
	float d = field(p).x; //min(min(sdTorus(rotMat1*vec3(c), vec2(10,3)),sdTorus(rotMat2*vec3(c), vec2(25,4))), -sdSphere(p, 50.0));
    
    #define CAPSULEDIST (0.0 + abs(sin(time))*1.0)
    //d = min(d, sdCapsule(rotMat1*vec3(c), vec3(0,CAPSULEDIST,0), vec3(0,-CAPSULEDIST,0), 4.0));
	return d < 0.0;
}

vec2 rotate2d(vec2 v, float a) {
	float sinA = sin(a);
	float cosA = cos(a);
	return vec2(v.x * cosA - v.y * sinA, v.y * cosA + v.x * sinA);	
}

#ifdef DITHERING
float dither(vec2 position, float brightness) {
	float bayer = texture(iChannel0, position).r;
    return step(bayer, brightness-0.1);
}
#endif

void main()
{
    vec2  screenPos    = squareFrame( resolution );
    float alpha = smoothstep( .2, 1., 1.- abs( screenPos.x ) );
    //uncomment below for fullscreen
    alpha = 1.;
    if( alpha <= 0. ) 
        discard;


	//vec2 uv = ((gl_FragCoord.xy * 2.0) / resolution.xy) - vec2(1);	// Make UV go from -1 to 1 instead of 0 to 1
    //uv.x *= resolution.x / resolution.y;

    /*
    vec3 s = vec3(sin(time*0.1)*45.0,sin(time*0.4)*15.0,cos(time*0.1)*45.0);
    #define FOCALLEN 0.6
    vec3 d = vec3(uv*FOCALLEN, 1.0);
    mat3 rotMat = rotationMatrix(vec3(0,1,sin(time*3.14159*0.1)*-0.3), -time*0.1 + 3.14159) * rotationMatrix(vec3(1,0,0), -0.4*sin(time*0.4) - 0.0);
    d = rotMat * d;
    
	vec3 rayDir = d;
	vec3 rayPos = s;
	*/

    mat3 rotMat = rotationMatrix(vec3(0,1,sin(time*3.14159*0.1)*-0.3), -time*0.1 + 3.14159) * rotationMatrix(vec3(1,0,0), -0.4*sin(time*0.4) - 0.0);

    vec3 rayDir = getRay( camera, target, screenPos, fov );
    vec3 rayPos = camera;
	
	ivec3 mapPos = ivec3(floor(rayPos + target * 0.1));

	vec3 deltaDist = abs(vec3(length(rayDir)) / rayDir) ;
	
	ivec3 rayStep = ivec3(sign(rayDir * vec3(.2)));

	vec3 sideDist = (sign(rayDir) * (vec3(mapPos) - rayPos) + (sign(rayDir) * 0.5) + 0.5) * deltaDist; 
	
	bvec3 mask = bvec3(0.);
    mat3 rotMat1 = rotationMatrix(vec3(1,1,0), 0.0/*time*0.3*/);
    mat3 rotMat2 = rotationMatrix(vec3(1,1,0), 0.0/*time*0.2*/);
	
	for (int i = 0; i < MAX_RAY_STEPS; i++) 
    {
		//if (getVoxel(mapPos)) continue;
		bvec3 b1 = lessThan(sideDist.xyz, sideDist.yzx);
		bvec3 b2 = lessThanEqual(sideDist.xyz, sideDist.zxy);
		mask.x = (b1.x && b2.x);
		mask.y = (b1.y && b2.y);
		mask.z = (b1.z && b2.z);

		//Would've done mask = b1 && b2 but the compiler is making me do it component wise.
		
		//All components of mask are false except for the corresponding largest component
		//of sideDist, which is the axis along which the ray should be incremented.			
		
        if(getVoxel(mapPos, rotMat1, rotMat2)) break;
		sideDist += vec3(mask)  * deltaDist;
		mapPos += ivec3(mask) * rayStep  ;
	}

    /*
		Basic lighting
		I calculate the distance from the current voxel center (mapPos) to a given light.
	*/
    
    gl_FragColor = vec4(0,0,0,1);	// Thanks otaviogood
    
    #define POW2(a) (a*a)
    
    #define CENTERCOLOR vec3(0,0.8,0.4)//(vec3(0,0.8,0.4) * clamp(cos(-time*2.0)*1.2-0.2, -0.1, 1.) )
    gl_FragColor.rgb += ( .05/POW2(distance(vec3(mapPos), rotMat*vec3(mapPos))) ) * 150.0 * CENTERCOLOR;
    
   /* #define MEDROTCOLOR vec3(0.1,0.5,0)
    rotMat = rotationMatrix(vec3(1,1,0), time*0.2);
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(sin(time)*25.0,0,cos(time)*25.0), rotMat*vec3(mapPos))) ) * 20.0 * MEDROTCOLOR;
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(sin(-time)*25.0,0,cos(time)*25.0), rotMat*vec3(mapPos))) ) * 20.0 * MEDROTCOLOR;
    
    #define CAPSULECOLOR (vec3(1,0,1)*(-cos(time*2.0)*0.5+0.5))
    //#define CAPSULEDIST (10.0 + sin(time)*5.0) Actually defined further up
    rotMat = rotationMatrix(vec3(1,1,0), time*0.3);
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(0, CAPSULEDIST+.5,0), rotMat*vec3(mapPos))) ) * 10.0 * CAPSULECOLOR;
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(0,-CAPSULEDIST+.5,0), rotMat*vec3(mapPos))) ) * 10.0 * CAPSULECOLOR;
    */
    //#define RIMCOLOR vec3(0,0.1,0.3) * max(0.0, sin(atan(float(mapPos.z), float(mapPos.x))*5.0+time*5.0)) * step(30.0, length(vec3(mapPos))) * (1.0-smoothstep(20., 50., abs(float(mapPos.y))))
    // /gl_FragColor.rgb += clamp(( 1.0/abs(sdTorus(vec3(mapPos - ivec3(0,0,0)), vec2(50.0,20)) )), 0., 1.0) * 5.0 * RIMCOLOR;
    
    /*#define OUTROTSPEED 0.2
    #define OUTROTRADIUS 45.0
    #define OUTROTBRIGHTNESS 100.0
    #define OUTROTCOLOR vec3(1,0.4,0)
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3( sin(time*OUTROTSPEED)*OUTROTRADIUS,0, cos(time*OUTROTSPEED)*OUTROTRADIUS), vec3(mapPos))) ) * OUTROTBRIGHTNESS * OUTROTCOLOR;
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3( cos(time*OUTROTSPEED)*OUTROTRADIUS,0,-sin(time*OUTROTSPEED)*OUTROTRADIUS), vec3(mapPos))) ) * OUTROTBRIGHTNESS * OUTROTCOLOR;
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(-sin(time*OUTROTSPEED)*OUTROTRADIUS,0,-cos(time*OUTROTSPEED)*OUTROTRADIUS), vec3(mapPos))) ) * OUTROTBRIGHTNESS * OUTROTCOLOR;
    gl_FragColor.rgb += ( 1.0/POW2(distance(vec3(-cos(time*OUTROTSPEED)*OUTROTRADIUS,0, sin(time*OUTROTSPEED)*OUTROTRADIUS), vec3(mapPos))) ) * OUTROTBRIGHTNESS * OUTROTCOLOR;
    */
    #ifdef DITHERING
    gl_FragColor.r = dither(fragCoord.xy / vec2(8), gl_FragColor.r);
    gl_FragColor.g = dither(fragCoord.xy / vec2(8), gl_FragColor.g);
    gl_FragColor.b = dither(fragCoord.xy / vec2(8), gl_FragColor.b);
    #endif

}