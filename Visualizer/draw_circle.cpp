#include <math.h>
#include <GL/glut.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
 
#define PI 3.1415926535897932384626433832795
 
#define WIDTH 1000
#define HEIGHT 1000
#define ITER_REPEAT 15
 
float pos[6] = {-1.0, -1.0, 1.0, 1.0, 30.0, -30.0};

void init();
void reshape(int, int);
void display();
 
//Specifies the radius of the circle
static float radius;
static int num_iters;
std::string buf;
std::ifstream sim_file;
int iteration;
std::vector< std::vector<float> > positions;

void preprocess() {
	int iter = 0;
	sim_file.open("simulation_data.txt");
	std::getline(sim_file, buf); 
	radius = std::stof(buf) * 100;
	std::getline(sim_file, buf); 
	num_iters = std::stoi(buf);
	positions.resize(num_iters);
	while (!sim_file.eof()) {
	    std::getline(sim_file, buf);
	    //std::cout << buf << "\n";
	    if (buf.length() == 0) { break; }
	    if (buf.compare("DONE WITH AN ITERATION") == 0) {
		//std::cout<<"LOL " << iter << "\n";
		iter++;
	    } else {
		positions[iter].push_back(std::stof(buf));
	    }
	}
	//std::cout  << "DONE \n";
	sim_file.close();
}

void update(int x) {
	if (iteration < num_iters) {
	    iteration++;
	    glutPostRedisplay();
	    glutTimerFunc(100, update, 0);
	} else {
	    std::cout << "DONE\n";
	}
}
 
int main(int argc, char **argv)
{
	preprocess();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutInitWindowPosition(50, 50);
	glutCreateWindow("OpenGL - Draw Circle");
	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutTimerFunc(100, update, 0);
	glutMainLoop();
 
	return 0;
}
 
void init()
{
	//Tells OpenGL to check which objects are in front of other objects
	//Otherwise OpenGL would draw the last object in front regardless
	//of it's position in Z space
	//Note: It is not necessary to enable this for a simple 2D circle
	//but is good practice
	
//glEnable(GL_DEPTH_TEST);
 
	//Tells OpenGL not to draw backfaces
	//Backfaces are defined by vertex drawing order
	//By default couter-clockwise drawing order specifies front faces
	//Note: The circle is drawn counter-clockwise
	//Note: It is not necessary to enable this for a simple 2D circle
	//but is good practice
//	glCullFace(GL_BACK);
//	glEnable(GL_CULL_FACE);
 
//	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
//	gluPerspective(45.0f, (GLfloat)WIDTH / (GLfloat)HEIGHT, 0.1f, 100000.0f);
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
 
	//Assign the clear screen color
	//Format (Red, Green, Blue, Alpha)
	//Values should remain normalized between 0 and 1
	glClearColor(1.0, 1.0, 1.0, 0.0);

}
 
void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)WIDTH, (GLsizei)HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(45.0f, (GLfloat)WIDTH / (GLfloat)HEIGHT, 0.1f, 100000.0f);
	glOrtho(0.0, 1000.0, 0.0, 1000.0, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
 
void display()
{
	//Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
	glPushMatrix();
        
 
		//Tells the camera where to be and where to look
		//Format (camera position x,y,z, focal point x,y,z, camera orientation x,y,z)
		//Remember that by default the camera points toward negative Z
	//	gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
 
	//	glPushMatrix();
 
			//Set Drawing Color - Will Remain this color until otherwise specified
			glColor3f(1.0, 0, 0);  //Some type of blue
 
			//Draw Circle
			//std::cout << iteration << " " << positions.size() << " " <<positions[iteration].size() << "\n";
			for(int j = 0; j < positions[iteration].size(); j+= 2) {
			glBegin(GL_POLYGON);
				//Change the 6 to 12 to increase the steps (number of drawn points) for a smoother circle
				//Note that anything above 24 will have little affect on the circles appearance
				//Play with the numbers till you find the result you are looking for
				//Value 1.5 - Draws Triangle
				//Value 2 - Draws Square
				//Value 3 - Draws Hexagon
				//Value 4 - Draws Octagon
				//Value 5 - Draws Decagon
				//Notice the correlation between the value and the number of sides
				//The number of sides is always twice the value given this range
			
				for(double i = 0; i < 2 * PI; i += PI / 12) //<-- Change this Value
 					glVertex3f(positions[iteration][j]*WIDTH + cos(i) * radius, positions[iteration][j+1]*HEIGHT+sin(i) * radius, 0.0);
					/*if (iteration > 0 && positions[iteration][j] != positions[iteration-1][j]) {
					    std::cout << "X POSITION IS DIFFERENT FOR " << j/2 << " AND THE DIFFERENCE IS " << (positions[iteration][j]-positions[iteration-1][j])*WIDTH<< "\n";
					}
					if (iteration > 0 && positions[iteration][j+1] != positions[iteration-1][j+1]) {
					    std::cout << "Y POSITION IS DIFFERENT FOR " << j/2 << " AND THE DIFFERENCE IS " << (positions[iteration][j+1]-positions[iteration-1][j+1])*WIDTH<< "\n";
					}*/
			glEnd(); }
			//std::cout << iteration;
			/*if (iteration < num_iters) {
			    iteration++;
			    num_times_per_iteration = 0;
			} */
			//Draw Circle
 
//		glPopMatrix();
 
	glPopMatrix();
 
	glFlush();
	glutSwapBuffers();
	//glFinish();
}
