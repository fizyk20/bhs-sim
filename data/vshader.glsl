varying vec4 position;

void main()
{
	position = gl_Vertex;
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
} 
