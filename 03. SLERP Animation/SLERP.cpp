#include <iostream>
#include <math.h>
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/include/GL/glut.h>

using namespace std;
using namespace Eigen;


// First window parameters
static int window_width = 800;
static int window_height = 600;
static int window_x = 100;
static int window_y = 50;
// Window dimension resize
static void on_reshape(int width,int height){
    window_height = height;
    window_width = width;
}

/*
// Camera position vars
static int cam_x = 10;
static int cam_y = 10;
static int cam_z = 5;
*/

// Euler's angles for the first and second state of the coordinate system
static double alpha_1 = M_PI/2, beta_1 = 3*M_PI / 4, gamma_1 = M_PI / 2;
static double alpha_2 = -M_PI/4, beta_2 = 2 * M_PI / 3, gamma_2 = -4*M_PI/3;
// Vertex of the coordinate system for the first and second state
static double x_1 = -1, y_1 = -2, z_1 = 6;
static double x_2 = 4, y_2 = 1, z_2 = -4;
// Time dependant parameters
static double alpha_t = 0, beta_t = 0, gamma_t = 0;
static double x_t = 0, y_t = 0, z_t = 0;
static Vector4d q_t;
// Quaternions of the first and last rotation of the system
static Vector4d q_1;
static Vector4d q_2;
// Time parameters
static double tm = 10;
static double t = 0;
static int timerActive = 0;
// Option parameter
static int ans;
// Matrix of rotation
static Matrix3d A;
double fi;


// Callback function stack
static void on_timer(int value);
static void on_keyboard(unsigned char key, int x, int y);
static void on_reshape(int width, int height);
static void on_display();

void Euler2A(double alpha, double beta, double gamma, Matrix3d &E2A){  
    // 
    // The rotation of the Euler's angles is given by the formula A =  
    
	Matrix3d RotXalpha;
	Matrix3d RotYbeta;
	Matrix3d RotZgamma;

	RotXalpha << 1,         0,           0, 
	             0, cos(alpha), -sin(alpha),
	             0, sin(alpha),  cos(alpha);

	RotYbeta << cos(beta), 0, sin(beta),
	            0,         1,         0,
	           -sin(beta), 0, cos(beta);

	RotZgamma << cos(gamma), -sin(gamma), 0,
	             sin(gamma),  cos(gamma), 0,
	             0,       0,              1;

	E2A = RotZgamma * RotYbeta * RotXalpha;
	return;
}
// AxisAngle[A] - returns eigenvector p = (px, py, pz) and the angle alpha(0, pi) so that the A = Rp(alpha)
void AxisAngle(Matrix3d &A, Vector3d &p, double &angle){
	// Is A*A.transpose() = E and is the determinant A.det = 1
	Matrix3d E;
	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	Matrix3d E_test;

	E_test = A.transpose() * A;

	Matrix3d E_test_2;
	E_test_2 << abs(round(E_test(0, 0))), abs(round(E_test(0, 1))), abs(round(E_test(0, 2))),
	            abs(round(E_test(1, 0))), abs(round(E_test(1, 1))), abs(round(E_test(1, 2))),
	            abs(round(E_test(2, 0))), abs(round(E_test(2, 1))), abs(round(E_test(2, 2)));

	if((E_test_2 != E) && (round(A.determinant()) != 1)){
		cout << "You have entered an incorrect matrix!" << endl;
		return;
	}

	// Evaluation of p - The vector around which we will rotate
    Vector3d v1;
	v1 << A(0, 0)-1, A(0, 1), A(0, 2);

	Vector3d v2;
	v2 << A(1, 0), A(1, 1)-1, A(1, 2);

	v1 = v1.normalized();
	v2 = v2.normalized();

	Vector3d zero_v;
	zero_v << 0, 0, 0;

	Vector3d crossp;
	crossp = v1.cross(v2);
	crossp(0) = round(crossp(0));
	crossp(1) = round(crossp(1));
	crossp(2) = round(crossp(2));

	// Are the vectors linearly dependant 
	if(crossp == zero_v){
		v2 << A(2, 0), A(2, 1), A(2, 2)-1;
		v2 = v2.normalized();
	}

	p = v1.cross(v2);
	p = p.normalized();

	// Finding the angle of rotation
	// By using vector p and a vector that is perpendicular to it (e.g. vectors v1 or v2)
	// we use matrix of rotation A which rotates vector v1 for the angle that we are calculating
	// and we use the formula that returns the angle between the two vectors

	Vector3d v1p;
	v1p = A * v1;
	
    angle = acos(v1.dot(v1p));

	// We check if we used the correct direction of p and if not we correct it
	Matrix3d mixp;
	mixp << p, v1, v1p;

	if(mixp.determinant() < 0){
		p = -1 * p;
	}

	return;
}
// Rodriguez[p, φ]  
// PRE: line p around which we rotate for the angle φ
// POST: Matrix of rotation around the oriented eigenvector p for the angle φ.
Matrix3d Rodriguez(Vector3d &p, double angle){
    // We get the matrix that we return from the following formula:
	// Rp(φ) = ppT + cos(φ)*(E − ppT) + sin(φ)*p×

	// Where p× is the matrix of the vector multiplication of the idenity vector's line p
	Matrix3d Px;
	Px <<  0,       -p(2, 0),  p(1, 0),
	       p(2, 0),  0,       -p(0, 0),
	      -p(1, 0),  p(0, 0),  0;

	// Identity matrix
	Matrix3d E;
	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	// Calculation of the Rodriguez matrix
	Matrix3d RodM;
	RodM = p*p.transpose() + cos(angle)*(E - p*p.transpose()) + sin(angle)*Px;

	return RodM;
}

// A2Euler[A]
// PRE: For the given orthogonally matrix A where det(A) = 1
// POST: Returns Euler's angles φ, θ, ψ, A = Rz(ψ) Ry(θ) Rx(φ)
void A2Euler(Matrix3d &RodM){
	double psi, teta, fi; 
		if (RodM(2, 0) < 1){
			if (RodM(2, 0) > -1){                       // Unique multiplication
				psi = atan2(RodM(1, 0), RodM(0, 0));
				teta = asin(-RodM(2, 0));
				fi = atan2(RodM(2, 1), RodM(2, 2));
			} else {                                    // Not unique, case: Ox3 = −Oz
				psi = atan2(-RodM(0, 1), RodM(1, 1));
				teta = M_PI/2;
				fi = 0;
			}
			} else {                                    // Not unique, case: Ox3 = Oz
				psi = atan2(-RodM(0, 1), RodM(1, 1));
				teta = -M_PI/2;
			    fi = 0;
			}

	return;
}


// AxisAngle2Q[p, φ]
// PRE: 
// POST: Returns unique quaternion q = (x, y, z, w) so that the Cq = Rp(φ) where p is a singular vector
void AxisAngle2Q(Vector3d& p, double angle, Vector4d& q){
	q << p(0, 0)*sin(angle/2), p(1, 0)*sin(angle/2), p(2, 0)*sin(angle/2), cos(angle/2);
	return;
}

// Q2AxisAngle[q]
// PRE: 
// POST: Returns a singular vector p = (px, py, pz) and the angle φ in range (o, π)
// so that the quaternion is a rotation Rp(φ), tj. Cq = Rp(φ) 
double QtoAxisAngle(Vector4d& q, Vector3d& vector_p){
	double angle_fi_temp = acos(q(3, 0));
	
    vector_p << q(0, 0)/sin(angle_fi_temp/2), q(1, 0)/sin(angle_fi_temp/2), q(2, 0)/sin(angle_fi_temp/2);

	angle_fi_temp = 2 * angle_fi_temp;

	vector_p = vector_p.normalized();

	return angle_fi_temp;
}


// SLERP algorithm
void SLERP(){

    // LERP in SLERP
    /*
    
    if (t < 0 || t > tm){
        cout << "T is not in valid range!" << endl;
        exit(0);
    }

    if (t == 0 || t == tm){
        q_t = q_1;
        return;
    }

    double q1_norma = q_1.norm();
    double q2_norma = q_2.norm();
    
    double cross = q_1.dot(q_2)/(q1_norma * q2_norma);
    if (cross < 0){  // Going on the shorter side of the sphere
        q_1 = -1*q_1;
        cross = -1*cross;
    }
    if (cross > 0.95) { // Quaternion is too close, bring back q_1
        q_t = q_1;
    }

    double fi = acos(cross);
    q_t = sin(fi*(1 - t/tm))/sin(fi)*q_1 + sin(fi*(t/tm))/sin(fi)*q_2;




    Vector4d zero_v;
	zero_v << 0, 0, 0, 0;
    Vector4d q_1n = q_1;
    Vector4d q_2n = q_2;

    Vector4d norm_q1 = q_1.normalized();
    double q1_norma = q_1.norm();
    if (norm_q1 != zero_v){
        q_1n(0) = q_1(0) / norm_q1(0);
        q_1n(1) = q_1(1) / norm_q1(1);
        q_1n(2) = q_1(2) / norm_q1(2);
        q_1n(3) = q_1(3) / norm_q1(3);
    }
    Vector4d norm_q2 = q_2.normalized();
    double q2_norma = q_2.norm();
    if (norm_q1 != zero_v){
        q_2n(0) = q_2(0) / norm_q2(0);
        q_2n(1) = q_2(1) / norm_q2(1);
        q_2n(2) = q_2(2) / norm_q2(2);
        q_2n(3) = q_2(3) / norm_q2(3);
    }

    double cos0 = q_1.dot(q_2);
    if (cos0 < 0){  // Going on the shorter side of the sphere
        q_1n = -q_1n;
        cos0 = -cos0;
    }
    if (cos0 > 0.95) { // Quaternion is too close, bring back q_1
        q_t = q_1n;
    }

    double fi0 = acos(cos0);
*/
    
    if(fi < M_PI/12 && fi > -M_PI/12){
        q_t = (1 - t/tm)*q_1 + (t/tm)*q_2;
        q_t = q_t.normalized();
    } else {
    //    q_t = ((sin(fi * (1 - t/tm))) / sin(fi)) * q_1 + ((sin(fi * (t/tm))) / sin(fi)) * q_2;
        q_t = sin(fi*(1 - t/tm))/sin(fi)*q_1 + sin(fi*(t/tm))/sin(fi)*q_2;
    }

    x_t = (1 - t/tm)*x_1 + (t/tm)*x_2;
    y_t = (1 - t/tm)*y_1 + (t/tm)*y_2;
    z_t = (1 - t/tm)*z_1 + (t/tm)*z_2;
    return;

}
// Coordinate system
void draw_coordinate_system(){

    GLUquadricObj *quadObj = gluNewQuadric();
	// (x, y, z)
    glPushMatrix();
    glColor3f(0, 1, 0);
    glTranslatef(3, 0, 0);
    glRotatef(90, 0, 1, 0);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0, 0);
    glTranslatef(0, 3, 0);
    glRotatef(-90, 1, 0, 0);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glPushMatrix();
    glColor3f(0, 0, 1);
    glTranslatef(0, 0, 3);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glColor3f(0.8, 0.8, 0.8);
    glutSolidSphere(0.1, 30, 30);

    glBegin(GL_LINES);
        glColor3f(0, 1, 0);
        glVertex3d(0, 0, 0);
        glVertex3d(3, 0, 0);
        
        glColor3f(1, 0, 0);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 3, 0);
        
        glColor3f(0, 0, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 0, 3);
    glEnd();
}

void draw_start_end(){
    Matrix3d E2A_1;

    Euler2A(alpha_1, beta_1, gamma_1, E2A_1);

    GLdouble matrix_1[16] = {E2A_1(0, 0), E2A_1(1, 0), E2A_1(2, 0), 0,
                          E2A_1(0, 1), E2A_1(1, 1), E2A_1(2, 1), 0,
                          E2A_1(0, 2), E2A_1(1, 2), E2A_1(2, 2), 0,
                          x_1, y_1, z_1, 1 };

    glPushMatrix();

    glMultMatrixd(matrix_1);
    draw_coordinate_system();

    glPopMatrix();

    Matrix3d E2A_2;

    Euler2A(alpha_2, beta_2, gamma_2, E2A_2);

    GLdouble matrix_2[16] = {E2A_2(0, 0), E2A_2(1, 0), E2A_2(2, 0), 0,
                          E2A_2(0, 1), E2A_2(1, 1), E2A_2(2, 1), 0,
                          E2A_2(0, 2), E2A_2(1, 2), E2A_2(2, 2), 0,
                          x_2, y_2, z_2, 1 };

    glPushMatrix();

    glMultMatrixd(matrix_2);
    draw_coordinate_system();

    glPopMatrix();

}

void calculate_quaternions(){
    // Calculating q1
    Matrix3d E2A_1;

    Euler2A(alpha_1, beta_1, gamma_1, E2A_1);

    Vector3d p_1;
    double angle_1;

    AxisAngle(E2A_1, p_1, angle_1);

    AxisAngle2Q(p_1, angle_1, q_1);

    // Calculating q2
    Matrix3d E2A_2;

    Euler2A(alpha_2, beta_2, gamma_2, E2A_2);

    Vector3d p_2;
    double angle_2;

    AxisAngle(E2A_2, p_2, angle_2);

    AxisAngle2Q(p_2, angle_2, q_2);

    return;
}

int main(int argc, char* argv[]){

    cout << "1 :: SLEPR (Spherical Linear Interpolation) algorithm" << endl;
    cout << "2 :: LERP (Linear Interpolation) algorithm" << endl;
    cin >> ans;
  
    // Input
    // 3 Euler's angles for the first state
    cout << "Input the 3 angles (alpha, beta, gamma) for the first state of the coordinate system:" << endl;
    cin >> alpha_1 >> beta_1 >> gamma_1;
    cout << "alpha: " << alpha_1 << " beta: " << beta_1 << " gamma: " << gamma_1 << endl;
    alpha_1 = alpha_1 * M_PI/180;
    beta_1 = beta_1 * M_PI/180;
    gamma_1 = gamma_1 * M_PI/180;

    // Center of the coordinate system for the first state
    cout << "Input the coordinates of the center of the coordinate system for the first state:" << endl;
    cin >> x_1 >> y_1 >> z_1;
    cout << "x: " << x_1 << " y: " << y_1 << " z: " << z_1 << endl;


    // 3 Euler's angles for the second state
    cout << "Input the 3 angles (alpha, beta, gamma) for the second state of the coordinate system:" << endl;
    cin >> alpha_2 >> beta_2 >> gamma_2;
    cout << "alpha: " << alpha_2 << " beta: " << beta_2 << " gamma: " << gamma_2 << endl;
    alpha_2 = alpha_2 * M_PI/180;
    beta_2 = beta_2 * M_PI/180;
    gamma_2 = gamma_2 * M_PI/180;

    // Center of the coordinate system for the second state
    cout << "Input the coordinates of the center of the coordinate system for the second state:" << endl;
    cin >> x_2 >> y_2 >> z_2;
    cout << "x: " << x_2 << " y: " << y_2 << " z: " << z_2 << endl;

    cout << "Enter the length of the animation: " << endl;
    cin >> tm;
    cout << "length: " << tm << endl;
    
    // Print simulation controls
    cout << endl << endl;
    cout << "--- Controls ---" << endl;
    cout << "S    ::: Start simulation" << endl;
    cout << "P    ::: Pause simulation" << endl;
    cout << "R    ::: Restart simulation" << endl;
    cout << "ESC  ::: Exit simulation" << endl;

    calculate_quaternions();

    double fi = acos(q_1(0)*q_2(0) + q_1(1)*q_2(1) + q_1(2)*q_2(2) + q_1(3)*q_2(3));
    if(fi > M_PI/2 || fi < -M_PI/2){
        q_1 = -q_1;
        fi = acos(q_1(0)*q_2(0) + q_1(1)*q_2(1) + q_1(2)*q_2(2) + q_1(3)*q_2(3)); 
    }
 
    // FreeGlut init
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

    // Window creation
    glutInitWindowPosition(window_x, window_y);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("SLERP animation");

    // Registering the callback functions
    glutKeyboardFunc(on_keyboard);
    glutReshapeFunc(on_reshape);
    glutDisplayFunc(on_display);

    // Setting the background and enabling depth test
    glClearColor(0, 0, 0, 0);
    glEnable(GL_DEPTH_TEST);

    // Line thickness
    glLineWidth(2);

    // Proram main loop entry
    glutMainLoop();

    return 0;
}

static void on_keyboard(unsigned char key, int x, int y){
    switch(key){
        case 27: // ESC program break
            exit(0);
            break;
        case 's':
        case 'S':
            // Start
            if (!timerActive) {
                glutTimerFunc(50, on_timer, 0);
                timerActive = 1;
            }
            break;
        case 'p':
        case 'P':
            // Pause
            timerActive = 0;
            break;
        case 'r':
        case 'R':            
            // Reset
            timerActive = 0;
            t = 0;
            break;
        /*
        case 37: // ArrowLeft
            cam_y -= 2;
            break;
        case 38: // ArrowUp
            cam_z += 2;
            break;
        case 39: // ArrowRight
            cam_y += 2;
            break;
        case 40: // ArrowDown
            cam_z -= 2;
            break;
        */
    }
}


static void on_timer(int value)
{
    // Callback == timer
    if (value != 0)
        return;

    // Simulation time update
    t += 0.1;

    // Simulation ending check
    if(t >= tm){
        t = 0;
        timerActive = 0;
        glutPostRedisplay();
        return;
    }

    // Redraw window
    glutPostRedisplay();

    // If needed timer reinitialized
    if (timerActive)
        glutTimerFunc(50, on_timer, 0);
}

// Play animation
static void on_display(){

    // Erasing previous buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Viewport adjusted
    glViewport(0, 0, window_width, window_height);

    // Projection adjusted
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, (float) window_width / window_height, 1, 200);

    // Camera position adjusted
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(15, 15, 5,
             0, 0, 0,
             0, 1, 0);

    // If needed matrix size:
    // GLfloat matrix[16]; 
    // glGetFloatv (GL_MODELVIEW_MATRIX , matrix);

    // Initial and last coordinate system
    draw_start_end();

    glPushMatrix();

    // SLEPR
    if(ans == 1){
        SLERP();
        
        // Moving quaternion to the rotation matrix
        Vector3d vector_p;
        double angle = QtoAxisAngle(q_t, vector_p);
    
        A = Rodriguez(vector_p, angle);
    }
    // LERP
    else if(ans == 2){
        alpha_t = (1 - t/tm)*alpha_1 + t/tm*alpha_2;
        beta_t = (1 - t/tm)*beta_1 + t/tm*beta_2;
        gamma_t = (1 - t/tm)*gamma_1 + t/tm*gamma_2;

        x_t = (1 - t/tm)*x_1 + t/tm*x_2;
        y_t = (1 - t/tm)*y_1 + t/tm*y_2;
        z_t = (1 - t/tm)*z_1 + t/tm*z_2;

        Euler2A(alpha_t, beta_t, gamma_t, A);
    }
    else {
        cout << "Not a valid option!" << endl;
        exit(1);
    }

    // Rotation and translation matrix
    GLdouble matrix[16] = {
                           A(0, 0), A(1, 0), A(2, 0), 0,
                           A(0, 1), A(1, 1), A(2, 1), 0,
                           A(0, 2), A(1, 2), A(2, 2), 0,
                               x_t,     y_t,     z_t, 1 
                          };

    glMultMatrixd(matrix);

    draw_coordinate_system();

    glPopMatrix();

    // Global coordinate system 
    glColor3f(0.8, 0.8, 0.8);
    glutSolidSphere(0.1, 30, 30);

    glBegin(GL_LINES);
        glColor3f(0, 1, 0);
        glVertex3d(-200, 0, 0);
        glVertex3d(200, 0, 0);
        
        glColor3f(1, 0, 0);
        glVertex3d(0, -200, 0);
        glVertex3d(0, 200, 0);
        
        glColor3f(0, 0, 1);
        glVertex3d(0, 0, -200);
        glVertex3d(0, 0, 200);
    glEnd();

    glutSwapBuffers();


}
