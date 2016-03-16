package ensembl

fun fromCSVString(CSVString: String, headers: Iterable<String>, separator: String = "\t"): Gene {
    val data = headers.zip(CSVString.split(separator)).map { it.first to it.second }.toMap()
    return Gene(
            id = data[ID_KEY] ?: "-1",
            start = data[START_KEY]?.toLong() ?: -1L,
            end = data[END_KEY]?.toLong() ?: -1L,
            strand = data[STRAND_KEY]?.toInt() ?: 0,
            chromosomeName = data[CHROMOSOME_NAME_KEY] ?: "-1"
    )
}

data class Gene(
        val id: String,
        val start: Long,
        val end: Long,
        val strand: Int,
        val chromosomeName: String
)


