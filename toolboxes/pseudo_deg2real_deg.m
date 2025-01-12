function [x_real, y_real]=pseudo_deg2real_deg(x_deg, y_deg, eye_x_deg, eye_y_deg)


d=42;
deg_pseudo=norm([x_deg, y_deg]);

x_cm=d*deg2rad(x_deg);
y_cm=d*deg2rad(y_deg);
eye_x_cm=d*deg2rad(eye_x_deg);
eye_y_cm=d*deg2rad(eye_y_deg);

OP=[x_cm+eye_x_cm, y_cm+eye_y_cm, d];
OQ=[eye_x_cm, eye_y_cm, d];

   OP_norm = OP / norm(OP);
   OQ_norm = OQ / norm(OQ);
   deg_real=rad2deg(acos(dot(OP_norm, OQ_norm))); 

   k=deg_real/deg_pseudo;

   x_real= k*x_deg;
   y_real= k*y_deg;

end