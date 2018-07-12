// Vertex shader for brick and mandelbrot shaders
// Derived from Orange Book Chapter 6

//  Light intensity and texture coords required by fragment shader
// varying float LightIntensity;
varying vec2 vTexCoord;
varying vec3 View;
varying vec3 Light;
varying vec3 Normal;

void main()
{
   /********** Pixel lighting **********/
   //  Vertex location in modelview coordinates
   vec4 P = gl_ModelViewMatrix * gl_Vertex;
   //  Light position
   Light  = gl_LightSource[0].position.xyz - P.xyz;
   //  Normal
   Normal = gl_NormalMatrix * gl_Normal;
   //  Eye position
   View  = -P.xyz;

   /********** Texture transparency **********/
   //  Get texture coords
   vTexCoord = gl_MultiTexCoord0.xy;

   //  Return fixed transform coordinates for this vertex
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
