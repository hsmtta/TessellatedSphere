# TessellatedSphere
The code generates triangulated sphere filled in almost equal area of triangle patches by tessellating icosahedron.

# Sample code

	void main()
	{
		const int nt = 10; // Tessellation number. Each icosahedron edge is devided by nt. 
		TessellatedSphere ts;
		ts.Init(nt);
		ts.SaveAsVtk(); // Save the result as VTK (ParaView) format.
	}
