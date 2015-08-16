#include "MyWorld.h"
#include "RigidBody.h"
#include "CollisionInterface.h"
#include "dart/dynamics/Skeleton.h"
#include "dart/dynamics/EllipsoidShape.h"
#include "dart/dynamics/BoxShape.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/Joint.h"
#include "dart/constraint/ConstraintSolver.h"
#include "dart/constraint/WeldJointConstraint.h"

using namespace Eigen;

MyWorld::MyWorld() {
  mFrame = 0;
  mTimeStep = 0.001;
  mGravity = Vector3d(0.0, -9.8, 0.0);
  mForce.setZero();
  // Create a collision detector
  mCollisionDetector = new CollisionInterface();
  epsilon = 1.0;
  // Create and intialize two default rigid bodies (You can add more rigid bodies if you want) 
  RigidBody *rb1 = new RigidBody(dart::dynamics::Shape::BOX, Vector3d(0.05, 0.05, 0.05));
  mCollisionDetector->addRigidBody(rb1, "box"); // Put rb1 in collision detector
  rb1->mPosition[0] = -0.3;
  rb1->mPosition[1] = -0.5;
  rb1->mAngMomentum = Vector3d(0.0, 0.1, 0.0);
  mRigidBodies.push_back(rb1);
    
  RigidBody *rb2 = new RigidBody(dart::dynamics::Shape::ELLIPSOID, Vector3d(0.06, 0.06, 0.06));
  mCollisionDetector->addRigidBody(rb2, "ellipse"); // Put rb2 in collision detector
  rb2->mPosition[0] = -0.1;
  rb2->mPosition[1] = -0.5;
  rb2->mAngMomentum = Vector3d(0.1, 0.0, 0.0);
  rb2->mColor = Vector4d(0.2, 0.8, 0.2, 1.0); // Blue
  mRigidBodies.push_back(rb2);

  RigidBody *rb3 = new RigidBody(dart::dynamics::Shape::ELLIPSOID, Vector3d(0.06, 0.06, 0.06));
  mCollisionDetector->addRigidBody(rb3, "ellipse"); // Put rb2 in collision detector
  rb3->mPosition[0] = -0.1;
  rb3->mPosition[1] = -0.3;
  rb3->mAngMomentum = Vector3d(0.1, 0.0, 0.0);
  rb3->mColor = Vector4d(0.2, 0.0, 0.6, 1.0); // Blue
  mRigidBodies.push_back(rb3);

  RigidBody *rb4 = new RigidBody(dart::dynamics::Shape::ELLIPSOID, Vector3d(0.06, 0.06, 0.06));
  mCollisionDetector->addRigidBody(rb4, "ellipse"); // Put rb2 in collision detector
  rb4->mPosition[0] = 0.3;
  rb4->mPosition[1] = -0.5;
  rb4->mAngMomentum = Vector3d(0.1, 0.0, 0.0);
  rb4->mColor = Vector4d(0.6, 0.1, 0.3, 1.0); // Blue
  mRigidBodies.push_back(rb4);

  RigidBody *rb5 = new RigidBody(dart::dynamics::Shape::BOX, Vector3d(0.05, 0.05, 0.05));
  mCollisionDetector->addRigidBody(rb5, "box"); // Put rb1 in collision detector
  rb5->mPosition[0] =   0.4;
  rb5->mPosition[1] = -0.5;
  rb5->mAngMomentum = Vector3d(0.0, 0.1, 0.0);
  rb5->mColor = Vector4d(0.9, 0.6, 0.0, 1.0); // Blue
  mRigidBodies.push_back(rb5);
}

void MyWorld::initializePinata() {
  // Add pinata to the collison detector
  mCollisionDetector->addSkeleton(mPinataWorld->getSkeleton(0));
  int nJoints = mPinataWorld->getSkeleton(0)->getNumBodyNodes();
  for (int i = 0; i < nJoints; i++) {
    int nDofs = mPinataWorld->getSkeleton(0)->getJoint(i)->getNumDofs();
    for (int j = 0; j < nDofs; j++)
      mPinataWorld->getSkeleton(0)->getJoint(i)->setDampingCoefficient(j, 1.0);
  }

  // Weld two seems to make a box
  dart::dynamics::BodyNode* top = mPinataWorld->getSkeleton(0)->getBodyNode("top");
  dart::dynamics::BodyNode* front = mPinataWorld->getSkeleton(0)->getBodyNode("front");
  dart::dynamics::BodyNode* back = mPinataWorld->getSkeleton(0)->getBodyNode("back");
  dart::constraint::WeldJointConstraint *joint1 = new dart::constraint::WeldJointConstraint(top, front);    
  dart::constraint::WeldJointConstraint *joint2 = new dart::constraint::WeldJointConstraint(top, back);    
  mPinataWorld->getConstraintSolver()->addConstraint(joint1);
  mPinataWorld->getConstraintSolver()->addConstraint(joint2);
}

MyWorld::~MyWorld() {
  for (int i = 0; i < mRigidBodies.size(); i++)
    delete mRigidBodies[i];
  mRigidBodies.clear();
  if (mCollisionDetector)
    delete mCollisionDetector;
  if (mPinataWorld)
    delete mPinataWorld;
}

void MyWorld::simulate() {
  mFrame++;
  // TODO: Replace the following code with your simulation
  for (int i = 0; i < mRigidBodies.size(); i++) {
        
        mRigidBodies[i] -> q.normalize();
       
        double w = mRigidBodies[i]->q.w();
        double x = mRigidBodies[i]->q.x();
        double y = mRigidBodies[i]->q.y();
        double z = mRigidBodies[i]->q.z();
        Matrix3d myMatrix;
        myMatrix.setIdentity();
        myMatrix << 1-2*pow(y,2)-2*pow(z,2), 2*x*y-2*w*z,     2*x*z+2*w*y,
                2*x*y+2*w*z,     1-2*pow(x,2)-2*pow(z,2), 2*y*z-2*w*x,
                2*x*z-2*w*y,     2*y*z+2*w*x,     1-2*pow(x,2)-2*pow(y,2);
        
        mRigidBodies[i] -> Iinv = myMatrix * mRigidBodies[i] -> Ibodyinv * (myMatrix.transpose());
        mRigidBodies[i] -> omega = mRigidBodies[i] -> Iinv * mRigidBodies[i] -> mAngMomentum;

        // mRigidBodies[i] -> Iinv = myMatrix * mRigidBodies[i] -> Ibody * (myMatrix.transpose());
        // mRigidBodies[i] -> omega = (mRigidBodies[i] -> Iinv).inverse() * mRigidBodies[i] -> mAngMomentum;
        Vector3d qtmp;
        qtmp.setZero();
        qtmp[0]=x; qtmp[1]=y; qtmp[2]=z;
        double first = 0.5 * (0 * w - mRigidBodies[i] -> omega.dot(qtmp));
        Vector3d qqdot = 0.5 * ( 0 * qtmp + w * mRigidBodies[i] -> omega + mRigidBodies[i] -> omega.cross( qtmp ));

        mRigidBodies[i] -> q.w() += first * mTimeStep;
        mRigidBodies[i] -> q.x() += qqdot[0] * mTimeStep; 
        mRigidBodies[i] -> q.y() += qqdot[1] * mTimeStep;
        mRigidBodies[i] -> q.z() += qqdot[2] * mTimeStep;

        // Quaterniond qtmp;
        // qtmp.setIdentity();
        // qtmp.x() = mRigidBodies[i] -> omega[0];
        // qtmp.y() = mRigidBodies[i] -> omega[1];
        // qtmp.z() = mRigidBodies[i] -> omega[2];

        // Quaterniond qqdot = 0.5 * qtmp * mRigidBodies[i]->q;

        // mRigidBodies[i]->q += qqdot * mTimeStep;
        //cout<<"quaternion "<<mRigidBodies[i] -> q.w()<<" "<<mRigidBodies[i] -> q.x()<<" "<<mRigidBodies[i] -> q.y()<<" "<<mRigidBodies[i] -> q.z()<<endl;
        w = mRigidBodies[i]->q.w();
        x = mRigidBodies[i]->q.x();
        y = mRigidBodies[i]->q.y();
        z = mRigidBodies[i]->q.z();

        mRigidBodies[i]->mOrientation << 1-2*pow(y,2)-2*pow(z,2), 2*x*y-2*w*z,     2*x*z+2*w*y,
                                        2*x*y+2*w*z,     1-2*pow(x,2)-2*pow(z,2), 2*y*z-2*w*x,
                                        2*x*z-2*w*y,     2*y*z+2*w*x,     1-2*pow(x,2)-2*pow(y,2);
        mRigidBodies[i]->mLinMomentum += mTimeStep * (mRigidBodies[i]->mMass * mGravity);
        mRigidBodies[i]->mVelocity = mRigidBodies[i]->mLinMomentum / mRigidBodies[i]->mMass;

        //mRigidBodies[i]->mPosition += mTimeStep * (mRigidBodies[i]->mLinMomentum / mRigidBodies[i]->mMass);
        
        mRigidBodies[i]->mPosition += mTimeStep * mRigidBodies[i]->mVelocity;
       
    }
  // Apply external force to the pinata
  mPinataWorld->getSkeleton(0)->getBodyNode("bottom")->addExtForce(mForce);
  mForce.setZero();
  // Simulate Pinata using DART
  mPinataWorld->step();
  // Run collision detector
  mCollisionDetector->checkCollision();

  // TODO: implement a collision handler
  collisionHandling();

  // Break the pinata if it has enough momentum
  if (mPinataWorld->getSkeleton(0)->getWorldCOMVelocity().norm() > 0.4)
    mPinataWorld->getConstraintSolver()->removeAllConstraints();
}

void MyWorld::collisionHandling() {
  int numContacts = mCollisionDetector->getNumContacts();
    for (int i = 0; i < numContacts; i++) 
    {
        RigidContact contact = mCollisionDetector->getContact(i);
        if ( contact.rb2 == NULL ) //Colliding with something static. Reflect
        { 
           
            Vector3d padot = (contact.rb1->mLinMomentum / contact.rb1->mMass) + (contact.rb1 -> omega.cross(contact.point - contact.rb1 -> mPosition ));
            Vector3d RA = contact.point - contact.rb1 -> mPosition;
            Matrix3d inersia;
            inersia << 0,0,0,
                        0,0,0,
                        0,0,0;
            Vector3d pbdot = contact.pinataVelocity;

            double vrel = (contact.normal).dot( padot - pbdot );

            
            if (vrel < 0.0 )
            {

              double numerator = -(1 + epsilon) * vrel;

                //calculate the denominator
              double term1 = 1 / contact.rb1 -> mMass;
              double term2 = 0;
              double term3 = contact.normal.dot((contact.rb1 -> Iinv * (RA.cross( contact.normal ))).cross(RA));
              double term4 = 0;

              double j = numerator / ( term1 + term2 + term3 + term4 );
              Vector3d J = j * contact.normal;

                //apply impluse to the bodies
              contact.rb1 -> mLinMomentum += J;
            
              contact.rb1 -> mAngMomentum += RA.cross(J);
              // contact.rb1 -> mVelocity = contact.rb1 -> mLinMomentum / contact.rb1 -> mMass;
              // contact.rb1 -> omega = contact.rb1 -> Iinv * contact.rb1 -> mAngMomentum;
            
            }
                
        }
        else
        {
            //Vector3d padot = contactPt_velocity( contact.rb1, contact.point);
            Vector3d padot = (contact.rb1->mLinMomentum / contact.rb1->mMass) + (contact.rb1 -> omega.cross(contact.point - contact.rb1 -> mPosition ));
            Vector3d pbdot = (contact.rb2->mLinMomentum / contact.rb2->mMass) + (contact.rb2 -> omega.cross(contact.point - contact.rb2 -> mPosition ));
            //Vector3d pbdot = contactPt_velocity( contact.rb2, contact.point);
            Vector3d RA = contact.point - contact.rb1 -> mPosition;
            Vector3d RB = contact.point - contact.rb2 -> mPosition;

            double vrel = (contact.normal).dot( padot - pbdot );

            if (vrel < 0.0 )
            {
              double numerator = -(1 + epsilon) * vrel;
            
            //calculate the denominator
              double term1 = 1 / contact.rb1 -> mMass;
              double term2 = 1 / contact.rb2 -> mMass;
              double term3 = contact.normal.dot((contact.rb1 -> Iinv * (RA.cross( contact.normal ))).cross(RA));
              double term4 = contact.normal.dot((contact.rb2 -> Iinv * (RB.cross( contact.normal ))).cross(RB));

              double j = numerator / ( term1 + term2 + term3 + term4 );
              Vector3d J = j * contact.normal;

            //apply impluse to the bodies
              contact.rb1 -> mLinMomentum += J;
              contact.rb2 -> mLinMomentum -= J;
              contact.rb1 -> mAngMomentum += RA.cross(J);
              contact.rb2 -> mAngMomentum -= RB.cross(J);

            //update v and omega
            //   contact.rb1 -> mVelocity = contact.rb1 -> mLinMomentum / contact.rb1 -> mMass;
            //   contact.rb2 -> mVelocity = contact.rb2 -> mLinMomentum / contact.rb2 -> mMass;

            // contact.rb1 -> omega = contact.rb1 -> Iinv * contact.rb1 -> mAngMomentum;
            // contact.rb2 -> omega = contact.rb2 -> Iinv * contact.rb2 -> mAngMomentum;
          }
        }
    }
}
