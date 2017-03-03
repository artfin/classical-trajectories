$$
J_x = J \cos \varphi \sin \theta \\
J_y = J \sin \varphi \sin \theta \\
J_z = J \cos \theta
$$


Для определения производных гамильтониана $\mathcal{H}(J_x, J_y, J_z)$ по углам, определяющим направление углового момента, используем цепное правило:


$$
\frac{\partial \mathcal{H}}{\partial \theta} = sum_{\alpha = x,y,z} \frac{\partial \mathcal{H}}{\partial J_\alpha} \frac{\partial J_\alpha}{\partial \theta} \\
\frac{\partial \mathcal{H}}{\partial \varphi} = sum_{\alpha = x,y,z} \frac{\partial \mathcal{H}}{\partial J_\alpha} \frac{\partial J_\alpha}{\partial \varphi} \\
\frac{\partial \mathcal{H}}{\partial J} = sum_{\alpha = x,y,z} \frac{\partial \mathcal{H}}{\partial J_\alpha} \frac{\partial J_\alpha}{\partial J} \\
$$

Подставив производные компонент углового момента по переменным $J, \varphi, \theta$, получаем следующую систему линейных уравнений, связывающие производные по старым и новым переменным.

$$
\begin{bmatrix}
\frac{\partial \mathcal{H}}{\partial \theta} \\
\frac{\partial \mathcal{H}}{\partial \varphi} \\
\frac{\partial \mathcal{H}}{\partial J}
\end{bmatrix}
=
\begin{bmatrix}
J \cos \varphi \cos \theta & J \sin \varphi \cos \theta & - J \sin \theta \\
- J \sin \varphi \sin \theta & J \cos \varphi \sin \theta & 0 \\
\cos \varphi \sin \theta & \sin \varphi \sin \theta & \cos \theta 
\end{bmatrix}
\begin{bmatrix}
\frac{\partial \mathcal{H}}{\partial J_x}
\frac{\partial \mathcal{H}}{\partial J_y}
\frac{\partial \mathcal{H}}{\partial J_z}
\end{bmatrix}
$$

