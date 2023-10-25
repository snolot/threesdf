uniform vec2 resolution;
uniform float time;
uniform float fov;
uniform float raymarchMaximumDistance;
uniform float raymarchPrecision;
uniform vec3 camera;
uniform vec3 target;

uniform samplerCube cubemap;
uniform vec3 anchors[15];

const vec3 backgroundColor = vec3(0.6, 0.7, 0.8);
const vec3 lightDir = vec3(-0.36, -0.48, 0.8);
const float shadow = 0.7;
const int maxIter = 100;

#define CAST_SHADOW

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

// Voxel Hit returns true if the voxel at pos should be filled.
bool voxelHit(vec3 pos) {
    
    // Fill in this function.
    vec3 hash = fract(pos * vec3(5.3983, 5.4427, 6.9371));
    hash += dot(hash, hash.yzx + 19.19);
    return length(pos) + 2.0 * fract((hash.x + hash.y) * hash.z) < 6.0 + 0.5 * sin(3.0 * time);

}

// Voxel Color returns the color at pos with normal vector norm.
vec3 voxelColor(vec3 pos, vec3 norm) {
    
    // Fill in this function.
    return mix(vec3(0.3, 0.5, 0.8), vec3(1.0, 0.4, 0.0), 0.5 * (length(floor(pos)) - 4.0));

}

////////////////////////////////////////////////////////////////////////////////
// Fill in the functions above.
// The engine below does not need to be modified.
////////////////////////////////////////////////////////////////////////////////

float castRay(vec3 eye, vec3 ray, out float dist, out vec3 norm) {
    vec3 pos = floor(eye);
    vec3 ri = 1.0 / ray;
    vec3 rs = sign(ray);
    vec3 ris = ri * rs;
    vec3 dis = (pos - eye + 0.5 + rs * 0.5) * ri;
    
    vec3 dim = vec3(0.0);
    for (int i = 0; i < maxIter; ++i) {
        if (voxelHit(pos)) {
            dist = dot(dis - ris, dim);
            norm = -dim * rs;
            return 1.0;
        }
    
        dim = step(dis, dis.yzx);
		dim *= (1.0 - dim.zxy);
        
        dis += dim * ris;
        pos += dim * rs;
    }

	return 0.0;
}

void main() {

    float dist;
    vec3 norm;
    float hit = castRay(camera, target, dist, norm);
    vec3 pos = camera + dist * target;

    vec3 color = voxelColor(pos - 0.001 * norm, norm);
    float shade = dot(norm, lightDir);
    
#ifdef CAST_SHADOW
    float illuminated = 1.0 - castRay(pos + 0.001 * norm, lightDir, dist, norm);
    float light = (1.0 + shadow * (illuminated * max(shade, 0.0) - 1.0)) * (1.0 - max(-shade, 0.0));
#else
    float light = (1.0 + shadow * (max(shade, 0.0) - 1.0)) * (1.0 - max(-shade, 0.0));    
#endif
    
    gl_FragColor = vec4(mix(backgroundColor, light * color, hit), 1.0);
}