	FILE* strmWrite = fopen("graph.dot", "w");
	write(strmWrite, g, DotDrawing());
	fclose(strmWrite);

	return 0;
}
