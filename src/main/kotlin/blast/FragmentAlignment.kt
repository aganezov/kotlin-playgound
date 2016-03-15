package blast


data class FragmentAlignment(
        val fragmentName: String = "",
        val subjectName: String = "chr-1",

        val fragmentLength: Long = -1L,

        val fragmentStrand: Int = -1,
        val subjectStrand: Int = -1,

        val fragmentStart: Long = -1L,
        val fragmentEnd: Long = -1L,

        val subjectStart: Long = -1L,
        val subjectEnd: Long = -1L,

        val gapsCnt: Long = -1L,
        val identity: Long = -1L,
        val alignmentLength: Long = -1L
        )