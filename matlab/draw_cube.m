function draw_cube(pos, size, color)
%{
 This functions draws a cube in 3D 
 pos : 3D coordinates of center of coube
 size: size/2 of side of the cube
 color: color of the cube
%}


%bottom rect
rec2_1x = [pos(1)-size, pos(1)+size, pos(1)+size, pos(1)-size];
rec2_1y = [pos(2)-size, pos(2)-size, pos(2)+size, pos(2)+size];
rec2_1z = [pos(3)-size, pos(3)-size, pos(3)-size, pos(3)-size];
patch(rec2_1x,rec2_1y,rec2_1z,color);
%top rect
rec2_2x = [pos(1)-size, pos(1)+size, pos(1)+size, pos(1)-size];
rec2_2y = [pos(2)-size, pos(2)-size, pos(2)+size, pos(2)+size];
rec2_2z = [pos(3)+size, pos(3)+size, pos(3)+size, pos(3)+size];
patch(rec2_2x,rec2_2y,rec2_2z, color);
%front rect 
rec2_3x = [pos(1)-size, pos(1)+size, pos(1)+size, pos(1)-size];
rec2_3y = [pos(2)-size, pos(2)-size, pos(2)-size, pos(2)-size];
rec2_3z = [pos(3)-size, pos(3)-size, pos(3)+size, pos(3)+size];
patch(rec2_3x,rec2_3y,rec2_3z,color);
%back rect
rec2_4x = [pos(1)-size, pos(1)+size, pos(1)+size, pos(1)-size];
rec2_4y = [pos(2)+size, pos(2)+size, pos(2)+size, pos(2)+size];
rec2_4z = [pos(3)-size, pos(3)-size, pos(3)+size, pos(3)+size];
patch(rec2_4x,rec2_4y,rec2_4z,color);
%right rect
rec2_5x = [pos(1)+size, pos(1)+size, pos(1)+size, pos(1)+size];
rec2_5y = [pos(2)-size, pos(2)+size, pos(2)+size, pos(2)-size];
rec2_5z = [pos(3)-size, pos(3)-size, pos(3)+size, pos(3)+size];
patch(rec2_5x,rec2_5y,rec2_5z,color);
%left rect
rec2_6x = [pos(1)-size, pos(1)-size, pos(1)-size, pos(1)-size];
rec2_6y = [pos(2)-size, pos(2)+size, pos(2)+size, pos(2)-size];
rec2_6z = [pos(3)-size, pos(3)-size, pos(3)+size, pos(3)+size];
patch(rec2_6x,rec2_6y,rec2_6z,color);

    






end