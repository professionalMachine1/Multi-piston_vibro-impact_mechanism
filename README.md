Работа посвящена исследованию динамики механизма с кривошипно-шатунным возбудителем колебаний при различном числе поршней-ударников.

В программе реализованы несколько основных методов, которые можно разделить на две части:
  - Методы, основанные на многократном решении дифференциального уравнения (MSV)
    * Осциллограмма (oscillogram_MSV)
    * Бифуркационная диаграмма (bifurcation_diagramm_MSV)
    * Карта динамических режимов (dynamic_mode_map_MSV)
  - Методы, основанные на точечном отображении (DMV)
    * Диаграмма Ламерея (point_map_DMV)
    * Бифуркационная диаграмма (bifurcation_diagramm_DMV)
    * Диаграмма Ляпунова (lyapunov_diagramm_DMV)
    * Карта динамических режимов (dynamic_mode_map_DMV)

Методы DMV работают в десятки раз быстрее методов MSV, но они ограничены в режимах работы. Такая разница в скорости обусловлена тем, что расчёты в методах DMV основаны 
на решении нелинейных уравнений, которые выполняются достаточно быстро. В свою очередь методы MSV основаны на решении дифференциального уравнения, которое просчитывается
довольно медленно, так как оно нелинейное, и шаг интегрирования довольно мал.

Чтобы строить диаграммы и карты для разных параметров используются макросы, которые определены в файле "SDE.h"

Литература:
  - Никифорова И. В. Динамика многопоршневых виброударных механизмов с кривошипно-шатунным возбудителем колебаний : дис. – Нижний Новгород : автореф. дис.… канд. физ.-мат. наук, 2017.
  - Килин А. А. Введение в теорию точечных отображений. Динамический хаос: учеб. пособие. – 2021.
  - Неймарк Ю. И. Метод точечных отображений в теории нелинейных колебаний. – URSS, 2010.
