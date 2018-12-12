radiation.dat
===============

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **timestep**
      - \-
      - time step counter
    * - **year**
      - \-
      - year
    * - **T**
      - \-
      - decimal day of year (DOY)
    * - **ShortIn (Rin)**
      - W m-2
      - Incoming shortwave radiation (copy from input)
    * - **LongIn (Rli)**
      - W m-2
      - Incoming longwave radiation (copy from input)
    * - **HemisOutShort (Eouto)**
      - W m-2
      - hemispherical outgoing shortwave radiation
    * - **HemisOutLong (Eoutt + Eoutte)**
      - W m-2
      - hemispherical outgoing longwave radiation
    * - **HemisOutTot (Eouto + Eoutt + Eoutte)**
      - W m-2
      - total hemispherical outgoing radiation
    * - **Net (Rntot)**
      - W m-2
      - total net radiation
