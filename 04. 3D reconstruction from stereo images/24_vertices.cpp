#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Geometry>
#include </usr/include/GL/glut.h>

using namespace std;
using namespace Eigen;

static int x_coord = 15;
static int y_coord = 45;
static int z_coord = 25;

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

MatrixXd reconstructed(24, 3);
        
MatrixXd calculate();
MatrixXd uAfine(MatrixXd &xx);
MatrixXd triD(MatrixXd &xx, MatrixXd &yy, MatrixXd &t1, MatrixXd &t2);

Vector3d cross_hidden(Vector3d & a, Vector3d & b, Vector3d & c, 
                             Vector3d & d, Vector3d & e, Vector3d & f, 
                             Vector3d & g, Vector3d & h, Vector3d & i, Vector3d & j);

void draw_small(MatrixXd &reconstructed);
void draw_medium(MatrixXd &reconstructed);
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
    gluPerspective(20,
                    window_width/(float)window_height, 1, 500);
    
}

static void on_display(void) {
    // Erasing previous buFFer
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
               0,  -1,  0
             );
    
    if(i == 0){
        reconstructed << calculate() * 0.5;
        cout << "Reconstructed vertices: " << endl << reconstructed << endl;
        i++;
    }
    glPushMatrix();
        draw_small(reconstructed);
        draw_medium(reconstructed);
        draw_big(reconstructed);
    glPopMatrix();
    glutSwapBuffers();
}


MatrixXd calculate() {
    MatrixXd x1(1, 3), x2(1, 3), x3(1, 3), x4(1, 3), x19(1, 3), x20(1, 3), x23(1, 3), x24(1, 3);
    MatrixXd y1(1, 3), y2(1, 3), y3(1, 3), y4(1, 3), y19(1, 3), y20(1, 3), y23(1, 3), y24(1, 3);
  
    // Eight vertices which will be used for the fundimental matrix F = FF
    
     x1 << 331,  75, 1;     
     y1 << 389,  76, 1;    
     x2 << 495,  53, 1;
     y2 << 561,  75, 1;
     x3 << 716, 166, 1;
     y3 << 566, 199, 1;
     x4 << 538, 191, 1;
     y4 << 371, 195, 1;
    x19 << 924, 600, 1;    
    y19 << 861, 655, 1;
    x20 << 700, 779, 1;    
    y20 << 456, 778, 1;
    x23 << 918, 787, 1;    
    y23 << 857, 838, 1;
    x24 << 697, 988, 1;    
    y24 << 461, 977, 1;
    

    MatrixXd xx(8, 3), yy(8, 3);
       
    xx<< x1, x2, x3, x4, x19, x20, x23, x24 ;
    yy<< y1, y2, y3, y4, y19, y20, y23, y24 ;

    // No normalization
    // y^T F x = 0
    // EQUATION [{a1, a2, a3}, {b1, b2, b3}]= { a1b1, a2b1, a3b1, a1b2, a2b2, a3b2, a1b3,a2b3,a3b3};
    // 8 equations given by the corespondencies 
    MatrixXd equation8 (8, 9);
    for (int i = 0 ;i < 8; i++) {
       int k = 0, l = 0;
        for (int j = 0 ; j < 9 ; j++) {
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

    MatrixXd SVDFFV  = SVDFF.matrixV();        
    MatrixXd SVDFFU  = SVDFF.matrixU();       
    MatrixXd SVDFFDD = SVDFF.singularValues();        
    
    // Third column is the missing epipole e1
    MatrixXd e1(1, 3);
    for(int i = 0; i < 3; i++) {
        e1(i) = SVDFFV(i, 2);
    }
    // Moving epipole e1 to afine
    for(int i = 0; i < 3; i++) {
        e1(i) = e1(i) / e1(2);
    }
     
    // Missing epipole2
    // e2 is the third column of U from the SVD
    MatrixXd e2(1, 3);
    for(int i = 0; i < 3; i++) {
        e2(i) = SVDFFU(i, 2);
    }
    // Moving epipole e2 to afine
    for(int i = 0; i < 3; i++) {
        e2(i) = e2(i) / e2(2);
    }
        
    // Making a diagonal matrix
    MatrixXd DD1 = SVDFFDD.array().matrix().asDiagonal();
    DD1(2,2) = 0;
    
    MatrixXd FF1 = SVDFFU * DD1 * (SVDFFV.transpose());    
    
    //
    // Reconstruct hidden vertices     
    //
    // Vectors of the hidden vertices:     

    Vector3d x_5, x_6, x_7, x_8, x_9, x_10, x_11, x_12, x_13, x_14, x_15, x_16, x_17, x_18, x_21, x_22;
    Vector3d y_5, y_6, y_7, y_8, y_9, y_10, y_11, y_12, y_13, y_14, y_15, y_16, y_17, y_18, y_21, y_22;
    
    // Making vectors out of the old vertices, for the cross product
    Vector3d x_1, x_2, x_3, x_4, x_19, x_20, x_23, x_24;
    Vector3d y_1, y_2, y_3, y_4, y_19, y_20, y_23, y_24;
    
    
    x_1  <<  x1(0),  x1(1),  x1(2);
    y_1  <<  y1(0),  y1(1),  y1(2);   
    x_2  <<  x2(0),  x2(1),  x2(2);
    y_2  <<  y2(0),  y2(1),  y2(2);       
    x_3  <<  x3(0),  x3(1),  x3(2);
    y_3  <<  y3(0),  y3(1),  y3(2);   
    x_4  <<  x4(0),  x4(1),  x4(2);
    y_4  <<  y4(0),  y4(1),  y4(2);       
    x_19 << x19(0), x19(1), x19(2);
    y_19 << y19(0), y19(1), y19(2);  
    x_20 << x20(0), x20(1), x20(2);
    y_20 << y20(0), y20(1), y20(2);  
    x_23 << x23(0), x23(1), x23(2);
    y_23 << y23(0), y23(1), y23(2);    
    x_24 << x24(0), x24(1), x24(2);
    y_24 << y24(0), y24(1), y24(2);
    
    x_5  << 330, 294, 1;
    x_7  << 714, 400, 1;
    y_7  << 567, 424, 1;
    x_8  << 553, 430, 1;
    y_8  << 378, 421, 1;   
    x_9  << 262, 339, 1;
    y_9  << 281, 311, 1;    
    y_10 << 713, 330, 1;
    x_11 << 774, 369, 1;
    y_11 << 688, 406, 1;
    x_12 << 312, 412, 1;
    y_12 << 234, 378, 1;
    x_13 << 262, 586, 1;
    y_14 << 718, 569, 1;
    x_15 << 770, 618, 1;
    y_15 << 686, 642, 1;
    x_16 << 312, 666, 1;    
    y_16 << 247, 615, 1;
    x_17 <<  91, 631, 1;    
    y_17 << 122, 552, 1;
    x_21 <<  95, 825, 1;    
    y_21 << 128, 721, 1;
 
    y_5 = cross_hidden(y_8, y_4, y_7, y_3, y_1,
                       y_4, y_1, y_3, y_2, y_8 );
    
    x_6 = cross_hidden(x_5, x_1, x_8, x_4, x_2,
                       x_8, x_5, x_3, x_2, x_7 );
    
    y_6 = cross_hidden(y_5, y_1, y_8, y_4, y_2, 
                       y_8, y_5, y_3, y_2, y_7 );
        

    x_10 = cross_hidden(x_16, x_13, x_12,  x_9, x_11,
                        x_12, x_11, x_16, x_15,  x_9 );
    
    y_13 = cross_hidden( y_15, y_16, y_10,  y_9, y_14,
                                y_16, y_12, y_15, y_11,  y_9 );
     
    x_14 = cross_hidden(x_16, x_15, x_12, x_11, x_13,
                        x_16, x_13, x_12,  x_9, x_15 ); 
    
    x_18 = cross_hidden(x_20, x_19, x_24, x_23, x_17, 
                        x_24, x_21, x_20, x_17, x_19 );
    
    y_18 = cross_hidden(y_20, y_19, y_24, y_23, y_17, 
                        y_24, y_21, y_20, y_17, y_19 );
    
    
    
    x_22 = cross_hidden(x_20, x_19, x_24, x_23, x_21,
                        x_24, x_21, x_20, x_17, x_23);
    
    y_22 = cross_hidden(y_20, y_19, y_24, y_23, y_21,
                        y_24, y_21, y_20, y_17, y_23); 
    
    
    MatrixXd x5(1, 3), x6(1, 3), x7(1, 3), x8(1, 3), x9(1, 3), x10(1, 3), x11(1, 3), x12(1, 3), x13(1, 3), x14(1, 3), x15(1, 3), x16(1, 3), x17(1, 3), x18(1, 3), x21(1, 3), x22(1, 3);
    MatrixXd y5(1, 3), y6(1, 3), y7(1, 3), y8(1, 3), y9(1, 3), y10(1, 3), y11(1, 3), y12(1, 3), y13(1, 3), y14(1, 3), y15(1, 3), y16(1, 3), y17(1, 3), y18(1, 3), y21(1, 3), y22(1, 3);
    
     x5 <<  x_5(0),  x_5(1),  x_5(2) ;
     x6 <<  x_6(0),  x_6(1),  x_6(2) ;
     x7 <<  x_7(0),  x_7(1),  x_7(2) ;
     x8 <<  x_8(0),  x_8(1),  x_8(2) ;    
     x9 <<  x_9(0),  x_9(1),  x_9(2) ;
    x10 << x_10(0), x_10(1), x_10(2) ;
    x11 << x_11(0), x_11(1), x_11(2) ;
    x12 << x_12(0), x_12(1), x_12(2) ;
    x13 << x_13(0), x_13(1), x_13(2) ;
    x14 << x_14(0), x_14(1), x_14(2) ;
    x15 << x_15(0), x_15(1), x_15(2) ;
    x16 << x_16(0), x_16(1), x_16(2) ;
    x17 << x_17(0), x_17(1), x_17(2) ;
    x18 << x_18(0), x_18(1), x_18(2) ;
    x21 << x_21(0), x_21(1), x_21(2) ;
    x22 << x_22(0), x_22(1), x_22(2) ;
    
    
    
     y5 <<  y_5(0),  y_5(1),  y_5(2) ; 
     y6 <<  y_6(0),  y_6(1),  y_6(2) ;
     y7 <<  y_7(0),  y_7(1),  y_7(2) ;
     y8 <<  y_8(0),  y_8(1),  y_8(2) ;
     y9 <<  y_9(0),  y_9(1),  y_9(2) ;
    y10 << y_10(0), y_10(1), y_10(2) ;
    y11 << y_11(0), y_11(1), y_11(2) ;
    y12 << y_12(0), y_12(1), y_12(2) ;
    y13 << y_13(0), y_13(1), y_13(2) ;
    y14 << y_14(0), y_14(1), y_14(2) ;
    y15 << y_15(0), y_15(1), y_15(2) ;
    y16 << y_16(0), y_16(1), y_16(2) ;
    y17 << y_17(0), y_17(1), y_17(2) ;
    y18 << y_18(0), y_18(1), y_18(2) ;
    y21 << y_21(0), y_21(1), y_21(2) ;
    y22 << y_22(0), y_22(1), y_22(2) ;
        
    
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
    MatrixXd picture1(24, 3), picture2(24, 3) ;
    picture1 << x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,
             x10, x11, x12, x13, x14, x15, x16, x17, x18,
             x19, x20, x21, x22, x23, x24 ;
    picture2 << y1,  y2,  y3,  y4,  y5,  y6,  y7,  y8,  y9,
             y10, y11, y12, y13, y14, y15, y16, y17, y18,
             y19, y20, y21, y22, y23, y24 ;
   
    MatrixXd equations(24, 4), tmp1(1, 3), tmp2(1, 3);
    for(int i = 0; i < 24; i++){
        tmp1 << picture1.row(i);
        tmp2 << picture2.row(i);
        equations.row(i) << triD(tmp1, tmp2, t1, t2);
    }
    // Translating coordinates to afine
    MatrixXd reconstructed(24, 3), pom(1, 4);
    for (int i = 0; i < 24; i++){
        pom << equations.row(i);
        reconstructed.row(i) << uAfine(pom); 
    }
    
    // We multiply the last coordinate e.g. 400, because it's so smll in comparison to the other coordinates, since 
    // we didin't normalize  
    MatrixXd reconstructed400(24, 3);
    for(int i = 0; i < 24; i++){
        reconstructed400.row(i) << reconstructed(i, 0), reconstructed(i, 1), reconstructed(i, 2)*400 ;
    }
 
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
    for(int i = 0; i < 4; i++){
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
void draw_small( MatrixXd &reconstructed) {
    MatrixXd small_edges (12, 2);
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
            glColor3f(0, 0, 1);
            glVertex3f(reconstructed(x, 0),reconstructed(x, 1), reconstructed(x, 2) );
            glVertex3f(reconstructed(y, 0),reconstructed(y, 1), reconstructed(y, 2));
        glEnd();
        
    }
    
}

// Medium box
void draw_medium(MatrixXd &reconstructed) {
    MatrixXd medium_edges (12, 2);
    medium_edges <<  8,  9,
                     9, 10,
                    10, 11,
                    11,  8,
                    12, 15,
                    15, 14,
                    14, 13,
                    13, 12,
                    15, 11,
                    14, 10,
                     8, 12,
                    13,  9;
    float x, y;
    for (int i = 0; i < 12; i++) {
        x = medium_edges(i, 0);
        y = medium_edges(i, 1);
        glBegin(GL_LINES);
            glColor3f(0, 1, 0);
            glVertex3f(reconstructed(x, 0),reconstructed(x, 1), reconstructed(x, 2) );
            glVertex3f(reconstructed(y, 0),reconstructed(y, 1), reconstructed(y, 2));
        glEnd();
        
    }
}

// Big box                
void draw_big( MatrixXd &reconstructed) {
    MatrixXd big_edges (12, 2);
    big_edges <<    16, 17, 
                    17, 18,
                    18, 19,
                    19, 16,
                    20, 21,
                    21, 22,
                    22, 23,
                    23, 20,
                    16, 20,
                    17, 21,
                    18, 22,
                    19, 23;
    
    float x, y;
    for (int i = 0; i < 12; i++) {
        x = big_edges(i, 0);
        y = big_edges(i, 1);
        glLineWidth(3.0f);
        glBegin(GL_LINES);
            glColor3f(1, 0, 0);
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
