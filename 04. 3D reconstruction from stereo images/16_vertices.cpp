#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Geometry>
#include </usr/include/GL/glut.h>

using namespace std;
using namespace Eigen;

static int x_coord = 30;
static int y_coord = 30;
static int z_coord = 15;

// First window parameters
static int window_width = 800;
static int window_height = 600;
static int window_x = 100;
static int window_y = 50;
// Window dimension resize

static int i = 0 ;

static void on_keyboard(unsigned char key, int x, int y);
static void on_reshape(int width, int height);
static void on_display(void);

MatrixXd reconstructed(16, 3);
        
MatrixXd calculate();
MatrixXd uAfine(MatrixXd &xx);
MatrixXd triD(MatrixXd &xx, MatrixXd &yy, MatrixXd &t1, MatrixXd &t2);

Vector3d cross_hidden(Vector3d & a, Vector3d & b, Vector3d & c, 
                             Vector3d & d, Vector3d & e, Vector3d & f, 
                             Vector3d & g, Vector3d & h, Vector3d & i, Vector3d & j);

void draw_small(MatrixXd &reconstructed);
void draw_big(MatrixXd &reconstructed);

void draw_axes();
int main (int argc, char *argv[]){
    cout << "Please input camera position coordinates (x, y, z):" << endl;
    cin >> x_coord >> y_coord >> z_coord;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("3D reconstruction");
	glClearColor(0, 0, 0, 0);
    
    glutKeyboardFunc(on_keyboard);
    glutReshapeFunc(on_reshape);
    glutDisplayFunc(on_display);

    glutMainLoop();

    return 0;
}


static void on_keyboard(unsigned char key, int x, int y){
    switch (key) {
        case 27:
            // Program exit
            exit(0);
            break;
    }
}

static void on_reshape(int width, int height) {
    window_width = width;
    window_height = height;
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40,
                    window_width/(float)window_height, 1, 500);
    
}

static void on_display(void) {
    // Erasing previous buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Viewport adjusted
    glViewport(0, 0, window_width, window_height);

    // Projection adjusted
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(20, (float) window_width / window_height, 1, 500);

    // Camera position adjusted
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(x_coord, y_coord, z_coord,
               0,  0,  0,
               0,  -1,  1
             );
    
    if(i == 0){
        reconstructed << calculate() * 0.05;
        cout << "Reconstructed vertices: " << endl << reconstructed << endl;
        i++;
    }
    glPushMatrix();
        draw_small(reconstructed);
        draw_big(reconstructed);
    glPopMatrix();
    glutSwapBuffers();
}



MatrixXd calculate() {
    MatrixXd x1(1, 3), x2(1, 3), x3(1, 3), x4(1, 3), x9(1, 3), x10(1, 3), x11(1, 3), x12(1, 3);
    MatrixXd y1(1, 3), y2(1, 3), y3(1, 3), y4(1, 3), y9(1, 3), y10(1, 3), y11(1, 3), y12(1, 3);
  
    // Eight vertices which will be used for the fundimental matrix F = FF

    x1  << 958,  38,   1;     
    y1  << 933,  33,   1;    
    x2  << 1117, 111,  1;
    y2  << 1027, 132,  1;
    x3  << 874,  285,  1;
    y3  << 692,  223,  1;
    x4  << 707,  218,  1;
    y4  << 595,  123,  1;
    x9  << 292,  569,  1;
    y9  << 272,  360,  1;    
    x10 << 770,  969,  1;
    y10 << 432,  814,  1;
    x11 << 770,  1465, 1;
    y11 << 414,  1284, 1;    
    x12 << 317,  1057, 1;    
    y12 << 258,  818,  1;
    

    MatrixXd xx(8, 3), yy(8, 3);
       
    xx<< x1, x2, x3, x4, x9, x10, x11, x12 ;
    yy<< y1, y2, y3, y4, y9, y10, y11, y12 ;
    
    // No normalization
    // y^T F x = 0
    // EQUATION [{a1, a2, a3}, {b1, b2, b3}]= { a1b1, a2b1, a3b1, a1b2, a2b2, a3b2, a1b3,a2b3,a3b3};
    // 8 equations given by the corespondencies 
    MatrixXd equation8 (8, 9);
    for (int i = 0; i < 8; i++) {
       int k = 0, l = 0;
        for (int j = 0; j < 9 ;j++) {
            equation8(i, j) = xx(i, k) * yy(i, l); 
            k++ ;
            if(k == 3){ 
                k = 0;
                l++;
            }   
        }
    }

    // DLT algorithm
    // SVD: 
    JacobiSVD<MatrixXd> SVDequation8(equation8, ComputeFullU | ComputeFullV);
    MatrixXd SVDequation8V = SVDequation8.matrixV();        
    MatrixXd p(1, 9);
    
    for(int i = 0; i < 9; i++){
        p(i) = SVDequation8V(i, 8);   
    }
       
    // Matrix form FF (3x3)
    MatrixXd FF(3, 3);
    FF << p(0), p(1), p(2),
          p(3), p(4), p(5),
          p(6), p(7), p(8);
        
    JacobiSVD<MatrixXd> SVDFF(FF, ComputeFullU | ComputeFullV);

    MatrixXd SVDffV  = SVDFF.matrixV();        
    MatrixXd SVDffU  = SVDFF.matrixU();       
    MatrixXd SVDffDD = SVDFF.singularValues();        
    
    // Third column is the missing epipole e1
    MatrixXd e1(1, 3);
    for(int i = 0; i < 3; i++) {
        e1(i) = SVDffV(i, 2);
    }
    // Moving epipole e1 to afine
    for(int i = 0; i < 3; i++) {
        e1(i) = e1 (i) / e1(2);
    }
     
    // Missing epipole2
    // e2 is the third column of U from the SVD
    MatrixXd e2(1, 3);
    for(int i = 0; i < 3; i++) {
        e2(i) = SVDffU(i, 2);
    }
    // Moving epipole e2 to afine
    for(int i = 0; i < 3; i++) {
        e2(i) = e2(i) / e2(2);
    }
        
    // Making a diagonal matrix
    MatrixXd DD1 = SVDffDD.array().matrix().asDiagonal();
    DD1(2,2) = 0;
    
    MatrixXd FF1 = SVDffU * DD1 * (SVDffV.transpose());    
    

    //
    // Reconstruct hidden vertices     
    //
    // Vectors of the hidden vertices:     
    Vector3d x_5, x_6, x_7, x_8, x_13, x_14, x_15, x_16;
    Vector3d y_5, y_6, y_7, y_8, y_13, y_14, y_15, y_16;
    
    // Making vectors out of the old vertices, for the cross product
    Vector3d x_1, x_2, x_3, x_4, x_9, x_10, x_11, x_12;
    Vector3d y_1, y_2, y_3, y_4, y_9, y_10, y_11, y_12;
    
    x_1  << x1(0),  x1(1), x1(2);
    y_1  << y1(0),  y1(1), y1(2);
    x_2  << x2(0),  x2(1), x2(2);
    y_2  << y2(0),  y2(1), y2(2);
    x_3  << x3(0),  x3(1), x3(2);
    y_3  << y3(0),  y3(1), y3(2);
    x_4  << x4(0),  x4(1), x4(2);
    y_4  << y4(0),  y4(1), y4(2);
    x_9  << x9(0),  x9(1), x9(2);
    y_9  << y9(0),  y9(1), y9(2);
    x_10 << x10(0), x10(1),x10(2);
    y_10 << y10(0), y10(1),y10(2);
    x_11 << x11(0), x11(1),x11(2);
    y_11 << y11(0), y11(1),y11(2);
    x_12 << x12(0), x12(1),x12(2);
    y_12 << y12(0), y12(1),y12(2);

    x_6 << 1094, 536, 1;
    y_6 << 980, 535, 1;
    x_7 << 862, 729, 1;
    y_7 << 652, 638, 1;
    x_8 << 710, 648, 1;
    y_8 << 567,532, 1;
    x_14 << 1487, 598, 1;
    y_14 << 1303, 700, 1;
    x_15 << 1462, 1079, 1;
    y_15 << 1257, 1165, 1;
    y_13 << 1077, 269, 1;
    
    x_5 = cross_hidden (x_4, x_8, x_6, x_2, x_1,
                        x_1, x_4, x_3, x_2, x_8 );
     
    y_5 = cross_hidden (y_4, y_8, y_6, y_2, y_1, 
                        y_1, y_4, y_3, y_2, y_8 );
        
    x_13 = cross_hidden(x_9, x_10, x_11, x_12, x_14,
                        x_11, x_15, x_10, x_14,  x_9 );
    
    
    x_16 = cross_hidden(x_10, x_14, x_11, x_15, x_12,
                        x_9, x_10, x_11, x_12, x_15 );
    
    y_16 = cross_hidden(y_10, y_14, y_11, y_15, y_12, 
                        y_9, y_10, y_11, y_12, y_15 );
    
    MatrixXd x5(1, 3), x6(1, 3), x7(1, 3), x8(1, 3), x13(1, 3), x14(1, 3), x15(1, 3), x16(1, 3);
    MatrixXd y5(1, 3), y6(1, 3), y7(1, 3), y8(1, 3), y13(1, 3), y14(1, 3), y15(1, 3), y16(1, 3);
    
    x5 <<  x_5(0),  x_5(1),  x_5(2)  ;
    x6 <<  x_6(0),  x_6(1),  x_6(2)  ;
    x7 <<  x_7(0),  x_7(1),  x_7(2)  ;
    x8 <<  x_8(0),  x_8(1),  x_8(2)  ;
    x13 << x_13(0), x_13(1), x_13(2) ;
    x14 << x_14(0), x_14(1), x_14(2) ;
    x15 << x_15(0), x_15(1), x_15(2) ;
    x16 << x_16(0), x_16(1), x_16(2) ;
    
    y5  << y_5(0),  y_5(1),  y_5(2)  ; 
    y6  << y_6(0),  y_6(1),  y_6(2)  ;
    y7  << y_7(0),  y_7(1),  y_7(2)  ;
    y8  << y_8(0),  y_8(1),  y_8(2)  ;
    y13 << y_13(0), y_13(1), y_13(2) ;
    y14 << y_14(0), y_14(1), y_14(2) ;
    y15 << y_15(0), y_15(1), y_15(2) ;
    y16 << y_16(0), y_16(1), y_16(2) ;
    
    
    
    //
    //  Triangulation
    //   

    MatrixXd m = MatrixXd::Identity(3, 3);
    MatrixXd t1 (3, 4);
    t1 << m.row(0), 0,
          m.row(1), 0,
          m.row(2), 0;  
          
    MatrixXd E2 (3, 3), E_2(3, 3);
    E2 <<     0, -e2(2),  e2(1),
          e2(2),      0, -e2(0),
         -e2(1),  e2(0),      0;
    
    E_2 = E2 * FF1;

    MatrixXd t2(3, 4);
    t2 << E_2.row(0), e2(0),
          E_2.row(1), e2(1),
          E_2.row(2), e2(2);  
          
    /* For every vertex we get a system of 4 equations with 4 homogenous missing variables.
    equations[x1, y1] = x1(2)*t1(3) -  x1(3)*t1(2), 
                       -x1(1)*t1(3) +  x1(3)*t1(1),
                        y1(2)*t2(3) -  y1(3)*t2(2),
                       -y1(1)*t2(3) +  y1(3)*t2(1);
    
    */
    
    // SVD for every vertex
    MatrixXd picture1(16, 3), picture2(16, 3) ;
    picture1 << x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
    picture2 << y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16 ;
   
    MatrixXd equations(16, 4), tmp1(1, 3), tmp2(1, 3);
    for(int i = 0; i < 16; i++){
        tmp1 << picture1.row(i);
        tmp2 << picture2.row(i);
        equations.row(i) << triD(tmp1, tmp2, t1, t2);
    }
    
    // Translating coordinates to afine
    MatrixXd reconstructed(16, 3), pom(1, 4);
    for (int i = 0; i < 16; i++){
        pom << equations.row(i);
        reconstructed.row(i) << uAfine(pom); 
    }
    
    // We multiply the last coordinate e.g. 400, because it's so smll in comparison to the other coordinates, since 
    // we didin't normalize
    MatrixXd reconstructed400(16, 3);
    
    reconstructed400 << reconstructed.col(0),reconstructed.col(1), reconstructed.col(2)*400;
     
    return reconstructed400;
}

MatrixXd uAfine(MatrixXd &xx) {
    MatrixXd afine(1, 3);
    afine << xx(0) / xx(3), xx(1) / xx(3), xx(2) / xx(3) ;
    
    return afine;
}


MatrixXd triD(MatrixXd &xx, MatrixXd &yy, MatrixXd &t1, MatrixXd &t2) {

    MatrixXd equation1(4, 4);
    equation1.row(0)<<      xx(1)*t1.row(2)  - xx(2)*t1.row(1);
    equation1.row(1)<< (-1)*xx(0)*t1.row(2)  + xx(2)*t1.row(0); 
    equation1.row(2)<<      yy(1)*t2.row(2)  - yy(2)*t2.row(1);
    equation1.row(3)<< (-1)*yy(0)*t2.row(2)  + yy(2)*t2.row(0);
    
    JacobiSVD<MatrixXd> svd_equation(equation1, ComputeFullU | ComputeFullV);

    MatrixXd svdV = svd_equation.matrixV();        
    MatrixXd svdU = svd_equation.matrixU();       
    MatrixXd svdDD = svd_equation.singularValues();        

    MatrixXd equation(1, 4);
    for( int i = 0 ; i < 4 ;i++){
        equation(i) = svdV(i, 3);
    }

    return equation;
}
    
void draw_axes(){
  glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0, 0, 0);
    glVertex3f(17, 0, 0);
  glEnd();
  glBegin(GL_LINES);
    glColor3f(0,1,0);
    glVertex3f(0,0,0);
    glVertex3f(0,17,0);
  glEnd();
  glBegin(GL_LINES);
    glColor3f(0,0,1);
    glVertex3f(0,0,0);
    glVertex3f(0,0,17);
  glEnd();
}

// Small box                 
void draw_small(MatrixXd &reconstructed) {
    
    MatrixXd small_edges (12, 2);
   
   /*
                 1, 2,
                 2, 3,
                 3, 4,
                 4, 1,
                 5, 6,
                 6, 7,
                 7, 8,
                 8, 5,
                 1, 5,
                 2, 6,
                 3, 7,
                 4, 8;
    */
   
    small_edges << 0, 1,
                 1, 2,
                 2, 3,
                 3, 0,
                 4, 5,
                 5, 6,
                 6, 7,
                 7, 4,
                 0, 4,
                 1, 5,
                 2, 6,
                 3, 7 ;
    
    float x, y;
    for (int i = 0; i < 12; i++) {
        x = small_edges(i, 0);
        y = small_edges(i, 1);
        glBegin(GL_LINES);
            glColor3f(1, 0, 0);
            glVertex3f(reconstructed(x, 0),reconstructed(x, 1), reconstructed(x, 2) );
            glVertex3f(reconstructed(y, 0),reconstructed(y, 1), reconstructed(y, 2));
        glEnd();
        
    }    
}

// Big box                
void draw_big(MatrixXd &reconstructed) {
    
    MatrixXd big_edges (12, 2);
    
    /*           
                  9, 10,
                 10, 11,
                 11, 12,
                 12, 9,
                 13, 14,
                 14, 15,
                 15, 16,
                 16, 13,
                  9, 13,
                  10, 14,
                 11, 15,
                 12, 16;
    */

    big_edges <<  8, 9,
                  9, 10,
                 10, 11,
                 11, 8,
                 12, 13,
                 13, 14,
                 14, 15,
                 15, 12,
                  8, 12,
                  9, 13,
                 10, 14,
                 11, 15;
    
    float x, y;
    for (int i = 0; i < 12; i++) {
        x = big_edges(i, 0);
        y = big_edges(i, 1);
        glBegin(GL_LINES);
            glColor3f(0, 1, 0);
            glVertex3f(reconstructed(x, 0),reconstructed(x, 1), reconstructed(x, 2) );
            glVertex3f(reconstructed(y, 0),reconstructed(y, 1), reconstructed(y, 2));
        glEnd();
        
    }    
}
Vector3d cross_hidden(Vector3d &a, Vector3d &b, Vector3d &c, 
                            Vector3d &d, Vector3d &e, Vector3d &f, 
                            Vector3d &g, Vector3d &h, Vector3d &i, Vector3d &j) {
    Vector3d result = (((a.cross(b)).cross(c.cross(d))).cross(e)).cross(((f.cross(g)).cross(h.cross(i))).cross(j)) ;
   
    result << result(0) / result(2), result(1) / result(2), result(2) / result(2);
    
    return result.array().round();
}
